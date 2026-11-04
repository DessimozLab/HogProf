"""Basic argument parsing, validation and builder dispatch checks."""

import io
import subprocess
import sys
import textwrap
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from types import ModuleType
from unittest.mock import Mock, patch

from hogprof import __version__, cli


class CliTest(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.work = Path(directory.name)
        self.tree = self.work / 'tree.nwk'
        self.tree.touch()
        self.orthoxml = self.work / 'family.orthoxml'
        self.orthoxml.touch()
        self.output = self.work / 'output'

        # Call the CLI without importing or executing the profile builder --
        # make a mock lshbuilder.run and substitute the real one
        builder_module = ModuleType('hogprof.lshbuilder')
        self.run = builder_module.run = Mock(return_value=('hashes', 'forest', 'matrix'))
        self.enterContext(patch.dict(sys.modules, {'hogprof.lshbuilder': builder_module}))
        self.startup = self.enterContext(patch.object(cli, 'print_startup'))
        self.enterContext(patch.object(cli, 'setup_cli'))

    def arguments(self):
        return [
            'build',
            '--orthoxml-glob', str(self.orthoxml),
            '--species-tree', str(self.tree),
            '--output-dir', str(self.output),
        ]

    def assert_cli_error(self, arguments, message):
        stderr = io.StringIO()
        with redirect_stderr(stderr), self.assertRaises(SystemExit) as error:
            cli.main(arguments)
        self.assertEqual(2, error.exception.code)
        self.assertIn(message, stderr.getvalue())
        self.run.assert_not_called()
        self.startup.assert_not_called()

    def test_missing_files_are_rejected(self):
        for source_flag in ('--oma', '--orthoxml-glob'):
            with self.subTest(source=source_flag):
                self.assert_cli_error([
                    'build',
                    source_flag, str(self.work / 'missing'),
                    '--tree', str(self.tree), '--output', str(self.output),
                ], source_flag)
        self.tree.unlink()
        self.assert_cli_error(self.arguments(), '--species-tree')

    def test_invalid_configurations(self):
        for options, message in (
            (['--jobs', '0'], '--jobs must be greater than zero'),
            (['--num-permutations', '0'], '--num-permutations must be greater than zero'),
            (['--min-species', '-1'], '--min-species must be nonnegative'),
            (['--min-events', '-1'], '--min-events must be nonnegative'),
            (['--loss-only', '--duplication-only'], 'not allowed with argument'),
            (['--reformat_names'], 'not supported'),
            (['--tax-weights', str(self.work / 'model')], '--tax-weights'),
        ):
            with self.subTest(options=options):
                self.assert_cli_error(self.arguments() + options, message)

    def test_output_file_is_rejected(self):
        self.output.touch()
        self.assert_cli_error(self.arguments(), '--output-dir: path is not a directory')

    def test_sources_are_mutually_exclusive(self):
        self.assert_cli_error([
            'build',
            '--tree', str(self.tree), '--output', str(self.output),
        ], 'one of the arguments')
        self.assert_cli_error(self.arguments() + ['--oma', str(self.orthoxml)],
                              'not allowed with argument')

    def test_build_dispatch(self):
        self.assertEqual(0, cli.main(self.arguments() + ['--jobs', '2', '--seed', '42']))
        self.run.assert_called_once_with(
            oma=None, orthoxml_glob=str(self.orthoxml), tree=self.tree,
            output_dir=self.output, tax_weights=None, taxon_mask=None,
            exclude_taxa=None, min_species=10, min_events=0, db_type=None,
            num_permutations=256, seed=42, loss_only=False,
            duplication_only=False, use_tax_ids=False, slice_subhogs=False,
            njobs=2, verbose=False,
        )
        self.startup.assert_called_once()

    def test_search_parameters_come_from_the_index(self):
        # --nsamples and --slice-subhogs are read from the manifest instead
        for option in (['--nsamples', '32'], ['--slice-subhogs']):
            with self.subTest(option=option):
                stderr = io.StringIO()
                with redirect_stderr(stderr), self.assertRaises(SystemExit) as error:
                    cli.main(['search', 'HOG:1', '--index-dir', str(self.work), *option])
                self.assertEqual(2, error.exception.code)
                self.assertIn('unrecognized arguments', stderr.getvalue())

    def test_process_arguments_are_used(self):
        with patch.object(sys, 'argv', ['hogprof'] + self.arguments()):
            self.assertEqual(0, cli.main())
        self.run.assert_called_once()

    def test_help(self):
        for arguments in ([], ['-h'], ['--help'], ['build', '--help'], ['search', '-h']):
            with self.subTest(arguments=arguments):
                stdout = io.StringIO()
                with redirect_stdout(stdout), self.assertRaises(SystemExit) as error:
                    cli.main(arguments)
                self.assertEqual(0, error.exception.code)
                self.assertIn('usage: hogprof', stdout.getvalue())
                if arguments and arguments[0] == 'build':
                    self.assertIn('usage: hogprof build', stdout.getvalue())
                    self.assertIn('--species-tree', stdout.getvalue())
                elif not arguments or arguments[0].startswith('-'):
                    self.assertIn('build', stdout.getvalue())
                    self.assertIn('search', stdout.getvalue())
        self.run.assert_not_called()

    def test_version(self):
        stdout = io.StringIO()
        with redirect_stdout(stdout), self.assertRaises(SystemExit) as error:
            cli.main(['--version'])
        self.assertEqual(0, error.exception.code)
        self.assertEqual(f'hogprof {__version__}\n', stdout.getvalue())
        self.run.assert_not_called()


class StartupImportsTest(unittest.TestCase):
    def test_help_and_version_avoid_execution_dependencies(self):
        # A fresh process catches eager imports hidden by the test runner's
        # module cache. Help must work without initializing the search stack.
        code = textwrap.dedent('''
            import sys
            from hogprof.cli import main

            try:
                main()
            except SystemExit as error:
                assert error.code == 0, error.code
            else:
                raise AssertionError("Expected help/version to exit")

            blocked = (
                "polars", "numpy", "pandas", "h5py", "tables", "datasketch",
                "pyoma", "pyham", "importlib.metadata", "pygments",
                "rich.progress", "rich.traceback", "rich.logging",
            )
            loaded = [name for name in sys.modules
                      if any(name == prefix or name.startswith(prefix + ".")
                             for prefix in blocked)]
            assert not loaded, loaded
            if sys.argv[1:] == ["--version"]:
                assert "rich" not in sys.modules
        ''')
        for arguments in ([], ['--help'], ['--version'], ['build', '--help'], ['search', '--help']):
            with self.subTest(arguments=arguments):
                completed = subprocess.run(
                    [sys.executable, '-B', '-c', code, *arguments],
                    capture_output=True, text=True,
                )
                self.assertEqual(0, completed.returncode, completed.stderr)


if __name__ == '__main__':
    unittest.main()
