"""User-facing console I/O for command-line tools"""

from __future__ import annotations
from abc import ABC, abstractmethod

import argparse
import sys
from collections.abc import Iterable
from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    import polars as pl
    from rich.console import Console

from hogprof import __version__
from hogprof.presets import DB_PRESETS

_Item = TypeVar("_Item")

# Rich Console shared by logging and progress, created only when needed
_console = None


def _get_console() -> Console:
    """Singleton for Console"""
    global _console
    if _console is None:
        from rich.console import Console
        _console = Console(stderr=True)
    return _console


def __getattr__(name):
    # Preserve access to the shared console without initializing Rich on import
    if name == "console":
        return _get_console()

    if name == "logger":
        import logging
        return logging.getLogger(__name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def setup_cli(verbose: bool = False,
              *,
              output_console: Console | None = None) -> None:
    """Configure application logging and terminal tracebacks for a CLI run."""
    import logging

    from rich.logging import RichHandler

    if output_console is None:
        output_console = _get_console()
    level = logging.DEBUG if verbose else logging.INFO

    handler = RichHandler(
        console=output_console,
        markup=False,
        rich_tracebacks=True,
        tracebacks_show_locals=verbose,
        show_path=verbose,
    )
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[handler],
        force=True,
    )

    # if in console (not redirected to file),
    # format error traceback nicely with Rich
    if output_console.is_terminal:
        from rich.traceback import install as install_rich_traceback
        install_rich_traceback(console=output_console,
                               show_locals=verbose,
                               word_wrap=True)

def print_startup(args=None, *, output_console: Console | None = None) -> None:
    """Print the banner and, when supplied, the initial CLI configuration."""
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text

    if output_console is None:
        output_console = _get_console()

    table = Table.grid(padding=(0, 3))
    table.add_column(style="dim")
    table.add_column()

    header = Table.grid(padding=(0, 1))
    header.add_column(justify="center")
    header.add_column()

    logo = Text("┌─┼─┐\n■ □ ■", style="bold cyan")
    name = Text.assemble(
        ("HogProf", "bold cyan"),
        (f" v{__version__}", "dim"),
    )

    header.add_row(logo, name)

    table.add_row(header)
    table.add_row()
    table.add_row(
        "Cite",
        "[link=https://doi.org/10.1371/journal.pcbi.1007553]"
        "Moi et al. (2020), PLOS Comp Biol[/link]",
    )
    if args is not None:
        source = (
            f"OMA: {args['oma']}" if args["oma"]
            else f"OrthoXML: {args['orthoxml_glob']}"
            if args["orthoxml_glob"]
            else f"OrthoXML tar: {args['orthoxml_tar']}"
        )
        table.add_section()
        table.add_row()
        table.add_row("Input", source)
        table.add_row("Output", str(args["output_dir"]))
        table.add_row("Workers", str(args["njobs"]))
        table.add_row("Seed", str(args["seed"]))

    output_console.print(
        Panel(table,
              title_align="left",
              subtitle_align="right",
              border_style="cyan",
              padding=(1, 2),
              expand=False,
        )
    )
    output_console.print("\n")
    if args is not None:
        output_console.print("[dim]Initializing...[/]")
        output_console.print("\n")


def track_progress(
    sequence: Iterable[_Item],
    *,
    description: str,
    total: int | None = None,
    output_console: Console | None = None,
) -> Iterable[_Item]:
    """Iterate with a Rich progress bar when output is an interactive terminal."""
    if output_console is None:
        output_console = _get_console()

    if output_console.is_terminal:
        from rich.progress import track
        # If output is CLI, it's all simple -- use Rich
        yield from track(
            sequence,
            description=description,
            total=total,
            console=output_console,
            disable=not output_console.is_terminal,
        )

    else:
        import logging
        import time
        from datetime import timedelta

        logger = logging.getLogger(__name__)
        # Otherwise if output is redirected to file (e.g. SLURM)
        # -- report progress differently
        if total is None:
            try:
                total = len(sequence)
            except TypeError:
                pass

        start = last_log = time.monotonic()
        last_percent = 0
        count = 0

        # how often in % to report
        log_step = 10
        # if log_step takes too long, how often in seconds to report anyway
        # -- every 1 min
        log_interval = 1 * 60

        first_message = f"{description}: 0% " + f"(0/{total})" if total else ""
        logger.info(first_message)
        for count, item in enumerate(sequence, 1):

            # First, yield item -- to feed the generator
            yield item

            now = time.monotonic()
            elapsed = now - start

            if total:
                percent = min(100, int(100 * count / total))
                due = percent >= last_percent + log_step
            else:
                percent = None
                due = True

            if due or now - last_log >= log_interval:
                if total:
                    eta = elapsed / count * (total - count)
                    eta_str = str(timedelta(seconds=round(eta)))
                    logger.info(
                        "%s: %d%% (%d/%d) | ETA %s",
                        description, percent, count, total, eta_str,
                    )
                    last_percent = percent
                else:
                    logger.info(
                        "%s: %d items | %.1f items/s",
                        description, count, count / elapsed,
                    )

                last_log = now

        logger.info(
            "%s: done (%d items in %.1fs)",
            description, count, time.monotonic() - start,
        )


def print_search_result(df: pl.DataFrame) -> None:
    """Pretty print for `hogprof search` results"""
    from rich.console import Console

    stdout = Console(file=sys.stdout)

    if not stdout.is_terminal:
        df.write_csv(sys.stdout, separator="\t")
        return

    from rich.table import Table

    table = Table(show_header=True, header_style="bold cyan", show_lines=False)

    for column in df.columns:
        table.add_column(
            column,
            justify="right" if df.schema[column].is_numeric() else "left",
        )

    for row in df.iter_rows():
        table.add_row(*(
            f"{value:.4f}" if isinstance(value, float)
            else "—" if value is None
            else str(value)
            for value in row
        ))

    stdout.print(table)


class RichArgumentParser(argparse.ArgumentParser):
    """Render argparse output with Rich while preserving streams and exit codes."""

    def print_help(self, file=None):
        from rich.console import Console

        file = sys.stdout if file is None else file
        print_startup(output_console=Console(file=file))
        super().print_help(file)

    def _print_message(self, message, file=None):
        if message:
            # Version output is plain text and needs no Rich rendering.
            if message == f"{self.prog} {__version__}\n":
                return super()._print_message(message, file)

            from rich.console import Console
            from rich.text import Text

            text = Text(message)
            text.highlight_regex(r"(?<!\w)--?[\w-]+", style="cyan")
            text.highlight_regex(
                r"(?m)^(?:usage|options|positional arguments):", style="bold"
            )
            if message.startswith(f"{self.prog}: error:"):
                text.stylize("bold red")
            Console(file=sys.stderr if file is None else file).print(
                text, end="", soft_wrap=True
            )

class ParserBuilder(ABC):
    """
    Arguments and validation for one CLI command.
    """
    def __init__(self, parser: argparse.ArgumentParser):
        self.parser = parser
        self.parsed_args: argparse.Namespace | None = None

    def parse(self, root_parser: argparse.ArgumentParser, argv: list[str]) -> argparse.Namespace:
        self.parsed_args = root_parser.parse_args(argv)
        return self.parsed_args

    @abstractmethod
    def validate(self) -> None:
        pass


class BuildParserBuilder(ParserBuilder):
    def __init__(self, parser: argparse.ArgumentParser):
        super().__init__(parser)
        self._add_args()

    def _add_args(self):
        from pathlib import Path
        parser = self.parser

        # Input files
        # - 'source' is one of three options: OMA DB, tarfile with .orthoxml, or glob
        source = parser.add_mutually_exclusive_group(required=True)
        source.add_argument(
            "--oma",
            "--OMA",  # compatibility option name; to remove
            dest="oma",
            help="Path to the OMA database file",
            type=Path,
        )
        source.add_argument(
            "--orthoxml-glob",
            "--OrthoGlob",  # compatibility option name; to remove
            dest="orthoxml_glob",
            type=str,
            help="A glob expression for orthoxml files",
        )
        source.add_argument(
            "--orthoxml-tar",
            "--tarfile",  # compatibility option name; to remove
            dest="orthoxml_tar",
            type=Path,
            help=argparse.SUPPRESS,
        )

        # - species tree
        parser.add_argument(
            "--species-tree",
            "--tree",
            "--mastertree",  # compatibility option name; to remove
            dest="tree",
            help="Newick species tree with node names matching the OrthoXML",
            type=Path,
            required=True,
        )

        # Output
        parser.add_argument(
            "--output-dir",
            "--output",
            "--outpath",  # compatibility option name; to remove
            "-o",
            dest="output_dir",
            help="Output directory path",
            type=Path,
            required=True,
        )

        # Configuration options
        parser.add_argument(
            "--tax-weights",
            "--taxweights",  # compatibility option name; to remove
            dest="tax_weights",
            help=argparse.SUPPRESS,
            type=str,
        )
        parser.add_argument(
            "--taxon-mask",
            "--taxmask",  # compatibility option name; to remove
            dest="taxon_mask",
            help="Keep only this clade, identified by its tree node name (e.g. Sauria)",
            type=str,
        )
        parser.add_argument(
            "--exclude-taxa",
            "--taxfilter",  # compatibility option name; to remove
            dest="exclude_taxa",
            help="Tree node names of clades to exclude",
            type=str,
            nargs="+",
        )
        parser.add_argument(
            "--min-species",
            "--specieslim",  # compatibility option name; to remove
            dest="min_species",
            help="Species threshold: root HOGs must exceed it; subHOGs must meet it (default:10)",
            type=int,
            default=10,
        )
        parser.add_argument(
            "--min-events",
            "--eventslim",  # compatibility option name; to remove
            dest="min_events",
            help="Minimum number of events (losses + duplications) a profile "
            "must have to be kept; 0 keeps all (default:0)",
            type=int,
            default=0,
        )
        parser.add_argument(
            "--db-type",
            "--dbtype",  # compatibility option name; to remove
            dest="db_type",
            help="Taxonomic range preset; --taxon-mask and --exclude-taxa override it",
            choices=DB_PRESETS,
        )
        parser.add_argument(
            "--num-permutations",
            "--nperm",  # compatibility option name; to remove
            dest="num_permutations",
            help="Number of hash functions used to construct each profile (default:256)",
            type=int,
            default=256,
        )
        parser.add_argument(
            "--seed",
            type=int,
            default=None,
            help="Random seed for reproducibility (default: random)",
        )

        # Flags
        events = parser.add_mutually_exclusive_group()
        events.add_argument(
            "--loss-only",
            "--lossonly",  # compatibility option name; to remove
            dest="loss_only",
            help="Only compile loss events",
            action="store_true",
        )
        events.add_argument(
            "--duplication-only",
            "--duplonly",  # compatibility option name; to remove
            dest="duplication_only",
            help="Only compile duplication events",
            action="store_true",
        )
        parser.add_argument(
            "--use-tax-ids",
            "--taxcodes",  # compatibility option name; to remove
            dest="use_tax_ids",
            help="Match species to tree nodes by NCBI taxon ID instead of species name",
            action="store_true",
        )
        parser.add_argument(
            "--slice-subhogs",
            "--slicesubhogs",  # compatibility option name; to remove
            dest="slice_subhogs",
            help="Make profiles for subHOGs",
            action="store_true",
        )

        # Retain unsupported legacy options only to give an actionable CLI error.
        parser.add_argument(
            "--reformat-names",
            "--reformat_names",
            dest="reformat_names",
            help=argparse.SUPPRESS,
            action="store_true",
        )

        # Other options
        parser.add_argument(
            "--jobs",
            "--njobs",
            "-j",
            "--nthreads",  # compatibility option name; to remove
            dest="njobs",
            help="Number of worker processes (default:1)",
            type=int,
            default=1,
        )

        parser.add_argument(
            "--verbose", "-v", help="Print verbose output", action="store_true"
        )
        parser.set_defaults(func=run_build)

    def validate(self):
        import glob
        from pathlib import Path
        parser = self.parser
        parsed_args = self.parsed_args

        if parsed_args.reformat_names:
            parser.error('--reformat-names is not supported in this version')

        if parsed_args.orthoxml_tar is not None:
            parser.error('--orthoxml-tar is not supported in this version; '
                         'extract the archive and use --orthoxml-glob')

        # check inputs
        for option, path in (('--oma', parsed_args.oma), ('--species-tree', parsed_args.tree)):
            if path is not None and not path.is_file():
                parser.error(f'{option}: input file does not exist or is not a file: {path}')

        # check Glob if provided
        if parsed_args.orthoxml_glob is not None:
            matches = glob.glob(parsed_args.orthoxml_glob)
            if not matches:
                parser.error('--orthoxml-glob: pattern does not match any input files')

            if any(not Path(path).is_file() for path in matches):
                parser.error('--orthoxml-glob: pattern must match only input files')


        if parsed_args.tax_weights is not None:
            for suffix in ('.json', '.h5'):
                path = Path(parsed_args.tax_weights + suffix)
                if not path.is_file():
                    parser.error(f'--tax-weights: model file does not exist or is not a file: {path}')

        if parsed_args.output_dir.exists() and not parsed_args.output_dir.is_dir():
            parser.error(f'--output-dir: path is not a directory: {parsed_args.output_dir}')

        for option, value in (('--jobs', parsed_args.njobs),
                              ('--num-permutations', parsed_args.num_permutations)):
            if value < 1:
                parser.error(f'{option} must be greater than zero')

        if parsed_args.min_species < 0:
            parser.error('--min-species must be nonnegative')

        if parsed_args.min_events < 0:
            parser.error('--min-events must be nonnegative')


class SearchParserBuilder(ParserBuilder):
    def __init__(self, parser: argparse.ArgumentParser):
        super().__init__(parser)
        self._add_args()

    def _add_args(self):
        parser = self.parser
        parser.add_argument("hog_id",
                            help="HOG ID for similarity search",
                            type=str)
        parser.add_argument("--k", "-k",
                            help="Number of nearest neighbors to consider (default:10)",
                            type=int, default=10)
        parser.add_argument("--level",
                            help="Taxonomic level at which the query is defined. "
                                 "Used for ambiguous HOG names defined at "
                                 "multiple levels.",
                            type=str)
        parser.add_argument("--index-dir",
                            help="Directory containing HogProf index files. "
                                 "Produced by `build`", type=str, required=True)
        parser.add_argument("--threshold", "-t",
                           help="Minimal significance level",
                           type=float,
                           default=None)
        parser.set_defaults(func=run_search)

    def validate(self):
        pass


class CalibrateParserBuilder(ParserBuilder):
    def __init__(self, parser: argparse.ArgumentParser):
        super().__init__(parser)
        self._add_args()

    def _add_args(self):
        from pathlib import Path
        from hogprof.presets import DEFAULT_ALPHAS

        parser = self.parser

        defaults = ", ".join(str(a) for a in DEFAULT_ALPHAS)
        parser.add_argument("--index-dir",
                            help="Directory containing HogProf index files. "
                                 "Produced by `build --slice-subhogs`",
                            type=Path, required=True)
        parser.add_argument("--threshold", "-t",
                            help=f"Significance levels to calibrate, in addition to {defaults}",
                            type=float, nargs="+", default=[])
        parser.add_argument("--sample-size",
                            help="Maximum number of subHOGs sampled per level (default:500)",
                            type=int, default=500)
        parser.add_argument("--seed",
                            type=int, default=None,
                            help="Random seed for reproducibility (default: read from the index)")
        parser.set_defaults(func=run_calibrate)

    def validate(self):
        parsed_args = self.parsed_args
        if parsed_args.sample_size < 2:
            self.parser.error("--sample-size must be at least 2")

        if any(not 0 < t < 1 for t in parsed_args.threshold):
            self.parser.error("--threshold values must be in (0, 1)")


def run_build(parsed_args):
    import logging

    setup_cli(verbose=parsed_args.verbose)

    # Suppress pyham's noisy INFO messages, including during initialization.
    logging.getLogger('pyham').setLevel(logging.WARNING)

    print_startup(vars(parsed_args))

    build_args = vars(parsed_args).copy()
    del build_args['reformat_names']
    del build_args['orthoxml_tar']
    del build_args['command']
    del build_args['func']

    # import and run as late as possible to not stagger CLI
    from hogprof.lshbuilder import run
    run(**build_args)


def run_search(parsed_args):
    setup_cli(verbose=False)

    search_args = vars(parsed_args).copy()
    del search_args["command"]
    del search_args["func"]

    from hogprof.profiler import run_search
    results = run_search(**search_args)
    print_search_result(results)


def run_calibrate(parsed_args):
    setup_cli()

    calibrate_args = vars(parsed_args).copy()
    del calibrate_args["command"]
    del calibrate_args["func"]

    from hogprof.background import run_calibrate
    run_calibrate(**calibrate_args)


def main(argv=None):
    # The parser renders its own help and errors. Logging and tracebacks are
    # configured only when dispatching an actual command.
    argv = sys.argv[1:] if argv is None else argv
    parser = RichArgumentParser(prog="hogprof")

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    # set up subparsers
    subparsers = parser.add_subparsers(dest="command", required=True)
    builder_parser = subparsers.add_parser("build", help="Build a searchable LSH index")
    search_parser = subparsers.add_parser("search", help="Query an LSH index")
    calibrate_parser = subparsers.add_parser("calibrate",
                                             help="Calibrate per-level significance thresholds for `search --threshold`")

    # Add arguments only for the requested command, to keep startup fast.
    # Anything else (no args, -h, --version, unknown command) is handled
    # by the top-level parser, which exits before validation.
    commands = {
        "build": (BuildParserBuilder, builder_parser),
        "search": (SearchParserBuilder, search_parser),
        "calibrate": (CalibrateParserBuilder, calibrate_parser),
    }
    if not argv:
        argv = ["-h"]
    elif argv == [argv[0]] and argv[0] in commands:
        # a bare command shows its help
        argv = [argv[0], "-h"]

    builder_class, command_parser = commands.get(argv[0], (None, None))
    if builder_class is None:
        parser.parse_args(argv)
        # unreachable: the subcommand is required, so argparse exits above
        return 2

    # parse, validate and run
    builder = builder_class(command_parser)
    args = builder.parse(parser, argv)
    builder.validate()
    args.func(args)
    return 0


if __name__ == '__main__':
    main()
