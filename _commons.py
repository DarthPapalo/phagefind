import logging
import sys
from pathlib import Path
from typing import Never

from rich.console import Console

PROGRAM_NAME: str = "PhageFind"
PROGRAM_DESC: str = (
    "A bioinformatic tool for the identification of bacteriophage sequences in de novo assemblies from bacteria sequencing data."
    # "\n\nPhageFind  Copyright (C) 2026  Pablo (DarthPapalo) Vidal\n"
    # "This program comes with ABSOLUTELY NO WARRANTY; for details see file LICENSE.\n"
    # "This is free software, and you are welcome to redistribute it\n"
    # "under certain conditions; see file LICENSE for details."
)
VERSION: str = "1.0.0"

# Global program rich console instance
console = Console()


def ensure_dir(
    logger: logging.Logger, path: Path, should_be_empty: bool, extra_msg: str = ""
) -> Path | Never:
    if path.exists():
        if not path.is_dir():
            logger.critical(
                f"'{path}' is not a valid directory{', ' + extra_msg if len(extra_msg) > 0 else ''}"
            )
            sys.exit(1)
        elif should_be_empty and len(list(path.iterdir())) > 0:
            logger.critical(
                f"'{path}' is not empty{', ' + extra_msg if len(extra_msg) > 0 else ''}"
            )
            sys.exit(1)

    path.mkdir(parents=True, exist_ok=True)

    return path
