import logging
import subprocess
from itertools import batched
from pathlib import Path
from typing import Iterator, cast

from program_config import PipelineConfig


def run_fastp(
    logger: logging.Logger,
    fastp_path: Path,
    config: PipelineConfig.PreprocessingConfig.PreprocessingFiltersConfig,
    reads_paths: list[Path],
    paired_ends: bool,
    results_path: Path,
) -> list[Path] | None:
    preprocessed_reads: list[Path] = []

    default_fastp_args: list[str] = [
        str(fastp_path),
        "--n_base_limit",
        str(config.n_base_limit),
        "--qualified_quality_phred",
        str(config.qualified_phred_quality),
        "--unqualified_percent_limit",
        str(config.unqualified_max_percentage),
        "--average_qual",
        str(config.minimum_average_quality),
        "--length_required",
        str(config.minimum_length),
        "--length_limit",
        str(config.maximum_length),
        "--report_title",
        "fastp report",
    ]

    # Create fastp reports folder, the program itself doesn't create it
    (results_path / "reports").mkdir()

    if paired_ends:
        for i, (read1, read2) in cast(
            Iterator[tuple[int, tuple[Path, Path]]],
            enumerate(batched(reads_paths, 2, strict=True)),
        ):
            read1_preprocessed: Path = results_path / f"preprocessed_{read1.name}"
            read2_preprocessed: Path = results_path / f"preprocessed_{read2.name}"
            fastp_paired_end_args: list[str] = [
                "--in1",
                str(read1),
                "--in2",
                str(read2),
                "--out1",
                str(read1_preprocessed),
                "--out2",
                str(read2_preprocessed),
                "--json",
                str(results_path / "reports" / f"qc_report_{i}.json"),
                "--html",
                str(results_path / "reports" / f"qc_report_{i}.html"),
            ]

            try:
                subprocess.run(
                    [*default_fastp_args, *fastp_paired_end_args],
                    stderr=subprocess.STDOUT,  # It prints normal output to stderr...
                    stdout=subprocess.DEVNULL,
                    check=True,
                )
            except subprocess.CalledProcessError:
                logger.critical(
                    f"fastp failed to preprocess paired end reads: {', '.join(map(str, reads_paths))}"
                )
                return None

            preprocessed_reads.extend((read1_preprocessed, read2_preprocessed))
    else:
        for i, read in cast(Iterator[tuple[int, Path]], enumerate(reads_paths)):
            preprocessed_read: Path = results_path / f"preprocessed_{read.name}"
            fastp_single_end_args: list[str] = [
                "-i",
                str(read),
                "-o",
                str(preprocessed_read),
                "--json",
                str(results_path / "reports" / f"qc_report_{i}.json"),
                "--html",
                str(results_path / "reports" / f"qc_report_{i}.html"),
            ]

            try:
                subprocess.run(
                    [*default_fastp_args, *fastp_single_end_args],
                    stderr=subprocess.STDOUT,
                    stdout=subprocess.DEVNULL,
                    check=True,
                )
            except subprocess.CalledProcessError:
                logger.critical(
                    f"fastp failed to preprocess reads: {', '.join(map(str, reads_paths))}"
                )
                return None

            preprocessed_reads.append(preprocessed_read)

    return preprocessed_reads
