"""User-facing console I/O for command-line tools"""

import logging
from collections.abc import Iterable, Iterator
from typing import TypeVar

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


def track_progress(
    sequence: Iterable[_Item],
    *,
    description: str,
    total: int | None = None,
    output_console: Console = console,
) -> Iterator[_Item]:
    """Iterate with a Rich progress bar when output is an interactive terminal."""
    return track(
        sequence,
        description=description,
        total=total,
        console=output_console,
        disable=not output_console.is_terminal)
