"""User-facing console I/O for command-line tools"""

import logging
from collections.abc import Iterable, Iterator
from typing import TypeVar
from HogProf import __version__

from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import track
from rich.traceback import install as install_rich_traceback

_Item = TypeVar("_Item")

# Rich Console shared by both logging and progress
console = Console(stderr=True)


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

def print_startup(args, *, output_console: Console = console) -> None:
    """Print the initial CLI configuration."""

    table = Table.grid(padding=(0, 3))
    table.add_column(style="dim")
    table.add_column()

    source = (
        f"OMA: {args['OMA']}" if args["OMA"]
        else f"OrthoXML: {args['OrthoGlob']}"
        if args["OrthoGlob"]
        else f"Tar: {args['tarfile']}"
    )

    table.add_row("Input", source)
    table.add_row("Output", str(args["outpath"]))
    table.add_row("Workers", str(args["njobs"]))

    title = Text.assemble((":: HogProf", "bold cyan"),
                          (f" v{__version__}", "dim"))

    output_console.print(
        Panel(table,
              title=title,
              subtitle="LSH index builder",
              title_align="left",
              subtitle_align="right",
              border_style="cyan",
              padding=(1, 2),
              expand=False,
        )
    )
    output_console.print("\n")
    output_console.print("[dim]Initializing...[/]")
    output_console.print("\n")


def track_progress(
    sequence: Iterable[_Item],
    *,
    description: str,
    total: int | None = None,
    output_console: Console = console,
) -> Iterator[_Item]:
    """Iterate with a Rich progress bar when output is an interactive terminal."""
    return track(sequence,
                 description=description,
                 total=total,
                 console=output_console,
                 disable=not output_console.is_terminal)
