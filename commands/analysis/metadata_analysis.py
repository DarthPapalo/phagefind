import csv
import logging
from collections import defaultdict
from pathlib import Path
from typing import Literal

from commands.analysis._commons import (
    METADATA_PHAGE_ACC_COLUMN,
    METADATA_PHAGE_HOST_COLUMN,
)


def metadata_analysis(
    logger: logging.Logger, phage_accs: set[str], metadata_files: list[Path]
) -> dict[str, set[str]]:
    """
    Returns a list of target organisms for each phage accession.
    """

    results: dict[str, set[str]] = defaultdict(set)
    for metadata_f in metadata_files:
        if not metadata_f.is_file():
            logger.error(
                f"Invalid metadata, '{str(metadata_f)}' is not a file, skipping"
            )
            continue
        if metadata_f.suffix.lower() not in [".tsv", ".csv"]:
            logger.error(
                f"Invalid metadata, '{str(metadata_f)}' is not a tsv or csv file, skipping"
            )
            continue

        delimiter: Literal["\t", ","] = (
            "\t" if metadata_f.suffix.lower() == ".tsv" else ","
        )
        metadata_reader: csv.DictReader[str] = csv.DictReader(
            metadata_f.read_text().splitlines(), delimiter=delimiter
        )
        for row in metadata_reader:
            if row[METADATA_PHAGE_ACC_COLUMN] in phage_accs:
                results[row[METADATA_PHAGE_ACC_COLUMN]].add(
                    row[METADATA_PHAGE_HOST_COLUMN]
                )

    return results
