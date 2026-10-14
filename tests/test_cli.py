"""Basic argument parsing, validation and builder dispatch checks."""

import io
import sys
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
            '--tree', str(self.tree), '--output', str(self.output),
        ], 'one of the arguments')
        self.assert_cli_error(self.arguments() + ['--oma', str(self.orthoxml)],
                              'not allowed with argument')


if __name__ == '__main__':
    unittest.main()
