import logging
from pathlib import Path
from typing import TypedDict

from program_config import PipelineConfig, ProgramsPaths

from .spades import run_spades


class AssemblyResults(TypedDict):
    """
    Results from the assembly step of the pipeline.
    """

    contigs_path: Path | None
    scaffolds_path: Path | None
    log_path: Path


def run(
    logger: logging.Logger,
    config: PipelineConfig.AssemblyConfig,
    programs_paths: ProgramsPaths,
    temp_path: Path,
    results_path: Path,
    reads_paths: list[Path],
    paired_ends: bool,
) -> AssemblyResults | None:
    """
    Generates a bacteria assembly from a set of genome reads.
    Returns results or log path in case of error.
    """
    if config.method == "SPAdes":
        workdir = temp_path / "spades"
        workdir.mkdir()

        spades_results: tuple[Path | None, Path | None, Path] | None = run_spades(
            logger, programs_paths.spades, config, reads_paths, paired_ends, workdir
        )

        contigs_path, scaffolds_path, log_path = spades_results

        return AssemblyResults(
            contigs_path=contigs_path.move_into(results_path)
            if contigs_path is not None
            else None,
            scaffolds_path=scaffolds_path.move_into(results_path)
            if scaffolds_path is not None
            else None,
            log_path=log_path.move_into(results_path),
        )
    else:
        logger.error(f"Unknown assembly method: {config.method}")
        return None
