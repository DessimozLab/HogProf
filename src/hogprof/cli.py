"""User-facing console I/O for command-line tools"""

import argparse
import glob
import logging
import sys
import time
from collections.abc import Iterable
from pathlib import Path
from typing import TypeVar

from rich.console import Console
from rich.logging import RichHandler
from rich.panel import Panel
from rich.progress import track
from rich.table import Table
from rich.text import Text
from rich.traceback import install as install_rich_traceback

from hogprof import __version__
from hogprof.presets import DB_PRESETS

_Item = TypeVar("_Item")

# Rich Console shared by both logging and progress
console = Console(stderr=True)
logger = logging.getLogger(__name__)


def setup_cli(verbose: bool = False,
              *,
              output_console: Console = console) -> None:
    """Configure application logging and terminal tracebacks for a CLI run."""
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
        install_rich_traceback(console=output_console,
                               show_locals=verbose,
                               word_wrap=True)

def print_startup(args=None, *, output_console: Console = console) -> None:
    """Print the banner and, when supplied, the initial CLI configuration."""

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
    output_console: Console = console,
) -> Iterable[_Item]:
    """Iterate with a Rich progress bar when output is an interactive terminal."""
    
    if output_console.is_terminal:
        # If output is CLI, it's all simple -- use Rich
        yield from track(
            sequence,
            description=description,
            total=total,
            console=output_console,
            disable=not output_console.is_terminal,
        )

    else:
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
                    logger.info(
                        "%s: %d%% (%d/%d) | ETA %.0fs",
                        description, percent, count, total, eta,
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



class RichArgumentParser(argparse.ArgumentParser):
    """Render argparse output with Rich while preserving streams and exit codes."""

    def print_help(self, file=None):
        file = sys.stdout if file is None else file
        print_startup(output_console=Console(file=file))
        super().print_help(file)

    def _print_message(self, message, file=None):
        if message:
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

def _validate_args(parser, parsed_args):
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

    #if parsed_args.min_events < 0:
    #    parser.error('--min-events must be nonnegative')


def main(argv=None):
    # setup Rich CLI before the parser. This is to provide
    # rich-formatted output for argparse argument validation, --help, etc.
    setup_cli()
    parser = RichArgumentParser(prog="lshbuilder")

    # Input files
    # - 'source' is one of three options: OMA DB, tarfile with .orthoxml, or glob
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--oma',
                        '--OMA', # compatibility option name; to remove
                        dest='oma',
                        help='Path to the OMA database file', type=Path)
    source.add_argument('--orthoxml-glob',
                        '--OrthoGlob', # compatibility option name; to remove
                        dest='orthoxml_glob',
                        type=str,
                        help='A glob expression for orthoxml files')
    source.add_argument('--orthoxml-tar',
                        '--tarfile', # compatibility option name; to remove
                        dest='orthoxml_tar',
                        type=Path,
                        help=argparse.SUPPRESS)

    # - species tree
    parser.add_argument('--species-tree', '--tree',
                        '--mastertree', # compatibility option name; to remove
                        dest='tree',
                        help='Newick species tree with node names matching the OrthoXML',
                        type=Path, required=True)

    # Output
    parser.add_argument('--output-dir', '--output',
                        '--outpath', # compatibility option name; to remove
                        '-o',
                        dest='output_dir',
                        help='Output directory path', type=Path, required=True)

    # Configuration options
    parser.add_argument('--tax-weights',
                        '--taxweights', # compatibility option name; to remove
                        dest='tax_weights',
                        help=argparse.SUPPRESS,
                        type=str
                        )
    parser.add_argument('--taxon-mask',
                        '--taxmask', # compatibility option name; to remove
                        dest='taxon_mask',
                        help='Keep only this clade, identified by its tree node name (e.g. Sauria)',
                        type=str)
    parser.add_argument('--exclude-taxa',
                        '--taxfilter', # compatibility option name; to remove
                        dest='exclude_taxa',
                        help='Tree node names of clades to exclude',
                        type=str, nargs='+')
    parser.add_argument('--min-species',
                        '--specieslim', # compatibility option name; to remove
                        dest='min_species',
                        help='Species threshold: root HOGs must exceed it; subHOGs must meet it',
                        type=int, default=10)
    parser.add_argument('--min-events',
                        '--eventslim', # compatibility option name; to remove
                        dest='min_events',
                        help='Keep subHOGs exceeding this loss or duplication count (-1 keeps all)',
                        type=int, default=0)
    parser.add_argument('--db-type',
                        '--dbtype', # compatibility option name; to remove
                        dest='db_type',
                        help='Taxonomic range preset; --taxon-mask and --exclude-taxa override it',
                        choices=DB_PRESETS)
    parser.add_argument('--num-permutations',
                        '--nperm', # compatibility option name; to remove
                        dest='num_permutations',
                        help='Number of hash functions used to construct each profile',
                        type=int, default=256)
    parser.add_argument('--seed',
                        type=int, default=None,
                        help="Random seed for reproducibility (default: random)")

    # Flags
    events = parser.add_mutually_exclusive_group()
    events.add_argument('--loss-only',
                        '--lossonly',  # compatibility option name; to remove
                        dest='loss_only',
                        help='Only compile loss events', action='store_true')
    events.add_argument('--duplication-only',
                        '--duplonly', # compatibility option name; to remove
                        dest='duplication_only',
                        help='Only compile duplication events', action='store_true')
    parser.add_argument('--use-tax-ids',
                        '--taxcodes', # compatibility option name; to remove
                        dest='use_tax_ids',
                        help='Match species to tree nodes by NCBI taxon ID instead of species name',
                        action='store_true')
    parser.add_argument('--slice-subhogs',
                        '--slicesubhogs', # compatibility option name; to remove
                        dest='slice_subhogs',
                        help='Make profiles for subHOGs', action='store_true')

    # Retain unsupported legacy options only to give an actionable CLI error.
    parser.add_argument('--reformat-names',
                        '--reformat_names',
                        dest='reformat_names',
                        help=argparse.SUPPRESS,
                        action='store_true')

    # Other options
    parser.add_argument('--jobs', '--njobs', '-j',
        '--nthreads', # compatibility option name; to remove
        dest="njobs",
        help="Number of worker processes", type=int, default=1)

    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')

    parser.add_argument('--verbose', '-v', help='Print verbose output', action='store_true')


    argv = sys.argv[1:] if argv is None else argv

    # print help if no args provided
    parsed_args = parser.parse_args(argv if argv else ['-h'])

    if parsed_args.verbose:
        setup_cli(verbose=True)

    _validate_args(parser, parsed_args)

    # Suppress pyham's noisy INFO messages, including during initialization.
    logging.getLogger('pyham').setLevel(logging.WARNING)

    print_startup(vars(parsed_args))

    build_args = vars(parsed_args).copy()
    del build_args['reformat_names']
    del build_args['orthoxml_tar']

    # import and run as late as possible to not stagger CLI
    from hogprof.lshbuilder import run
    run(**build_args)
    return 0


if __name__ == '__main__':
    main()
