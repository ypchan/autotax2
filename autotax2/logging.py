from __future__ import annotations

import sys
from contextlib import contextmanager
from typing import Iterator, Optional

try:
    from rich.console import Console
    from rich.logging import RichHandler
    from rich.panel import Panel
    from rich.progress import (
        BarColumn,
        Progress,
        SpinnerColumn,
        TaskProgressColumn,
        TextColumn,
        TimeElapsedColumn,
        TimeRemainingColumn,
    )
    from rich.table import Table
    RICH_AVAILABLE = True
except Exception:  # pragma: no cover - fallback for minimal environments
    Console = None  # type: ignore
    Progress = None  # type: ignore
    Panel = None  # type: ignore
    Table = None  # type: ignore
    RICH_AVAILABLE = False


class PlainConsole:
    def print(self, *args, **kwargs) -> None:
        print(*args)

    def rule(self, title: str = "") -> None:
        print(f"\n--- {title} ---")


console = Console(stderr=True) if RICH_AVAILABLE else PlainConsole()


def info(message: str) -> None:
    console.print(f"[cyan]INFO[/cyan] {message}" if RICH_AVAILABLE else f"INFO {message}")


def success(message: str) -> None:
    console.print(f"[green]OK[/green] {message}" if RICH_AVAILABLE else f"OK {message}")


def warn(message: str) -> None:
    console.print(f"[yellow]WARN[/yellow] {message}" if RICH_AVAILABLE else f"WARN {message}")


def error(message: str) -> None:
    console.print(f"[red]ERROR[/red] {message}" if RICH_AVAILABLE else f"ERROR {message}")


def print_kv(title: str, rows: dict[str, str]) -> None:
    if RICH_AVAILABLE:
        table = Table(title=title, show_lines=False)
        table.add_column("Key", style="bold")
        table.add_column("Value")
        for k, v in rows.items():
            table.add_row(str(k), str(v))
        console.print(table)
    else:
        print(f"== {title} ==")
        for k, v in rows.items():
            print(f"{k}\t{v}")


def print_example(title: str, command: str, notes: Optional[str] = None) -> None:
    if RICH_AVAILABLE:
        body = f"[bold]Command[/bold]\n\n[green]{command}[/green]"
        if notes:
            body += f"\n\n[bold]Notes[/bold]\n{notes}"
        console.print(Panel(body, title=title, border_style="cyan"))
    else:
        print(f"\n{title}\n{command}\n")
        if notes:
            print(notes)


@contextmanager
def progress_context() -> Iterator[object]:
    """Return a Rich Progress instance or a no-op fallback."""
    if RICH_AVAILABLE:
        progress = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=console,
            transient=False,
        )
        with progress:
            yield progress
    else:
        class NoProgress:
            def add_task(self, description, total=1):
                print(description)
                return 0

            def update(self, task_id, advance=1, completed=None, description=None):
                if description:
                    print(description)
        yield NoProgress()


def progress_step(progress: object, task_id: int, message: str, advance: int = 1) -> None:
    if RICH_AVAILABLE and hasattr(progress, "update"):
        progress.update(task_id, advance=advance, description=message)
    else:
        print(message)
