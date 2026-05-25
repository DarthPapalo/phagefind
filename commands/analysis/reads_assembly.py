import logging
from pathlib import Path
from typing import cast

from _commons import console, ensure_dir
from program_config import PipelineConfig, ProgramsPaths

from .modules.assembly import AssemblyResults
from .modules.assembly import run as run_assembly
from .modules.preprocessing import run as run_preprocessing


def reads_assembly(
    logger: logging.Logger,
    temp: Path,
    loaded_config: PipelineConfig,
    loaded_programs_paths: ProgramsPaths,
    query_path: Path,
    output_path: Path,
    reads_paths: list[Path],
    preprocessing: bool,
    paired_ends: bool,
) -> Path | None:
    """
    Returns the Assembly file path
    """

    # Check input reads
    for read_path in reads_paths:
        if not read_path.is_file():
            logger.error(f"Reads file '{read_path}' not found")
            return None

    logger.info(f"Starting reads analysis for: {', '.join(map(str, reads_paths))}")

    # ====== Preprocessing ======
    ready_reads_paths: list[Path] = []
    if preprocessing:
        preprocessing_path: Path = ensure_dir(
            logger,
            output_path / loaded_config.preprocessing.sub_dir,
            True,
            "can't run preprocessing of reads",
        )

        with console.status("Preprocessing reads...", spinner="dots"):
            ready_reads_paths: list[Path] | None = run_preprocessing(
                logger,
                loaded_config.preprocessing,
                loaded_programs_paths,
                preprocessing_path,
                reads_paths,
                paired_ends,
            )

            if ready_reads_paths is None:
                logger.error("Couldn't preprocess reads")
                return None

        logger.info("Reads preprocessing completed")

    else:
        ready_reads_paths = reads_paths

    # ====== Assembly ======
    assembly_results_path: Path = ensure_dir(
        logger,
        output_path / loaded_config.assembly.sub_dir,
        True,
        "can't run assembly step",
    )

    assembly_path: Path
    with console.status("Assembling reads...", spinner="dots"):
        assembly_results: AssemblyResults | None = run_assembly(
            logger,
            loaded_config.assembly,
            loaded_programs_paths,
            temp,
            assembly_results_path,
            ready_reads_paths,
            paired_ends,
        )
        if assembly_results is None:
            return None
        elif (
            assembly_results["contigs_path"] is None
            or assembly_results["scaffolds_path"] is None
        ):
            logger.error(
                f"Assembly couldn't complete, check '{assembly_results['log_path']}'"
            )
            return None

        assembly_path = cast(
            Path, assembly_results[loaded_config.assembly.use_file + "_path"]
        )

    logger.info("Reads assembly completed")

    return assembly_path
