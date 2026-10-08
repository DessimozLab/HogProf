"""User-facing console I/O for command-line tools"""

import argparse
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

from HogProf import __version__

DB_PRESETS = {
    'all': {'taxfilter': None, 'taxmask': None},
    'plants': {'taxfilter': None, 'taxmask': 33090},
    'archaea': {'taxfilter': None, 'taxmask': 2157},
    'bacteria': {'taxfilter': None, 'taxmask': 2},
    'eukarya': {'taxfilter': None, 'taxmask': 'Eukaryota'},
    'protists': {'taxfilter': [2, 2157, 33090, 4751, 33208], 'taxmask': None},
    'fungi': {'taxfilter': None, 'taxmask': 4751},
    'metazoa': {'taxfilter': None, 'taxmask': 33208},
    'vertebrates': {'taxfilter': None, 'taxmask': 7742},
}

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
            f"OMA: {args['OMA']}" if args["OMA"]
            else f"OrthoXML: {args['OrthoGlob']}"
            if args["OrthoGlob"]
            else f"Tar: {args['tarfile']}"
        )
        table.add_section()
        table.add_row()
        table.add_row("Input", source)
        table.add_row("Output", str(args["outpath"]))

        njobs = args["nthreads"] or args["njobs"]
        table.add_row("Workers", str(njobs))

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


def main(argv=None):
    # setup Rich CLI before the parser. This is to provide
    # rich-formatted output for argparse argument validation, --help, etc.
    setup_cli()
    parser = RichArgumentParser(prog="lshbuilder")

    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')
    parser.add_argument('--taxweights', help='load optimised weights from keras model',type = str)
    parser.add_argument('--taxmask', help='consider only one branch (e.g. Sauria)',type = str)
    parser.add_argument('--taxfilter', help='remove these taxa', type = str, nargs='*')
    parser.add_argument('--outpath', '-o', help='Output directory path', type=Path, required=True)
    parser.add_argument('--dbtype', help='preconfigured taxonomic ranges', choices=DB_PRESETS)
    parser.add_argument('--OMA', help='use oma data ', type = str)
    parser.add_argument('--OrthoGlob', help='a glob expression for orthoxml files ' , type = str)
    parser.add_argument('--tarfile', help='use tarfile with orthoxml data', type = str)
    parser.add_argument('--nperm', help='number of hash functions to use when constructing profiles',
                        type=int, default=256)
    parser.add_argument('--mastertree', help='master taxonomic tree. nodes should correspond to orthoxml' , type = str)

    # limits
    parser.add_argument('--specieslim', help='minimum number of species in a subhog' , type = int, default=10)
    parser.add_argument('--eventslim', help='minimum number of events (loss/duplication) in a subhog' , type = int, default=0)

    # multiprocessing
    parser.add_argument('--nthreads', help='[deprecated] Number of threads for multiprocessing', type=int)
    parser.add_argument("--njobs", help="Number of jobs for multiprocessing", type=int, default=1)

    # Flags
    parser.add_argument('--lossonly', help='only compile loss events', action='store_true')
    parser.add_argument('--duplonly', help='only compile duplication events', action='store_true')
    parser.add_argument('--taxcodes', help='use taxid info in HOGs', action='store_true')
    parser.add_argument('--reformat_names',
                        help='Correct broken species trees by replacing all names with numbers.',
                        action='store_true')
    parser.add_argument('--slicesubhogs', help='Make profiles for subhogs', action='store_true')
    parser.add_argument('--verbose', '-v', help='print verbose output', action='store_true')

    argv = sys.argv[1:] if argv is None else argv

    # print help if no args provided
    parsed_args = parser.parse_args(argv if argv else ['-h'])

    if parsed_args.verbose:
        setup_cli(verbose=True)

    if not (parsed_args.OMA or parsed_args.OrthoGlob or parsed_args.tarfile):
        parser.error('Please specify input data with --OMA, --OrthoGlob or --tarfile')

    if parsed_args.reformat_names:
        parser.error('--reformat_names is not supported in this version')

    # Suppress pyham's noisy INFO messages, including during initialization.
    logging.getLogger('pyham').setLevel(logging.WARNING)

    print_startup(vars(parsed_args))

    # import and run as late as possible to not stagger CLI
    from HogProf.lshbuilder import run
    return run(parsed_args)


if __name__ == '__main__':
    main()
