import logging
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import cast

import clapy

from _commons import console, ensure_dir
from commands.analysis._commons import (
    METADATA_PHAGE_ACC_COLUMN,
    METADATA_PHAGE_HOST_COLUMN,
    AnalysisMode,
    AnalysisResults,
)
from commands.analysis.metadata_analysis import metadata_analysis
from commands.analysis.modules.differences_analysis import (
    DifferencesAnalysisResults,
    HitsPerFeature,
)
from commands.analysis.modules.differences_analysis import (
    run as run_differences_analysis,
)
from commands.analysis.modules.identification import IdentificationResults
from commands.analysis.modules.identification import run as run_identification
from commands.analysis.reads_assembly import reads_assembly
from commands.analysis.report_generation import generate_report
from program_config import PipelineConfig, ProgramsPaths

ANALYSIS_COMMAND_NAME = "analyze"

# fmt: off
analysis_common_args: list[clapy.Arg] = [
    clapy.Arg("query-path")
        .help(
            "Multi-FASTA file containing the reference sequences to search in the bacteria assembly"
        )
        .long("--query")
        .value_parser(Path),
    # UNIMPLEMENTED ONLINE MODE
    # clapy.Arg("online-mode")
    #     .help("UNIMPLEMENTED - Enables BLAST searches for additional bacteria targets for the identified bacteriophages")
    #     .long("--online-mode")
    #     .flag(),
    clapy.Arg("metadata")
        .help(f"One or several metadata files for multiple phage target identification (TSV or CSV format and must include the '{METADATA_PHAGE_ACC_COLUMN}' and '{METADATA_PHAGE_HOST_COLUMN}' columns).")
        .long("--metadata")
        .nargs("*")
        .value_parser(Path)
        .default(())
        .append(),
]
# fmt: on

cli: clapy.Command = clapy.Command(ANALYSIS_COMMAND_NAME).help(
    "Perform a full analysis from sets of bacteria genome reads or an assembly."
)

# fmt: off
reads_analysis_command: clapy.Command = (
    clapy.Command("reads")
        .help("Perform a full analysis from a set of bacteria genome reads.")
        .arguments([
                clapy.Arg("reads-paths")
                    .help("List of bacteria genome reads.")
                    .nargs("+")
                    .value_parser(Path),
                clapy.Arg("preprocessing")
                    .help("Enables the preprocessing of the reads.")
                    .long("--preprocessing")
                    .flag(),
                clapy.Arg("paired-ends")
                    .help("Use if each pair of reads files corresponds to paired ends.")
                    .long("--paired-ends")
                    .flag(),
        ])
)
# fmt: on
cli.subcommand(reads_analysis_command.arguments(analysis_common_args))

# fmt: off
asembly_analysis_command: clapy.Command = (
    clapy.Command("assembly")
        .help("Perform the full analysis from a bacteria assembly.")
        .arguments([
                clapy.Arg("assembly-path")
                    .help("Bacterial assembly file path.")
                    .value_parser(Path),
        ])
)
# fmt: on
cli.subcommand(asembly_analysis_command.arguments(analysis_common_args))


def check_paths(
    logger: logging.Logger,
    output_path: Path,
):
    # Check output dir
    output_path.mkdir(exist_ok=True)
    if not output_path.is_dir():
        logger.critical(
            f"Output dir '{output_path}' is not a valid directory, aborting pipeline."
        )
        sys.exit(1)
    elif len(list(output_path.iterdir())) > 0:
        logger.critical(
            f"Output directory '{output_path}' is not empty, aborting pipeline."
        )
        sys.exit(1)


def run_analysis(
    mode: AnalysisMode,
    logger: logging.Logger,
    config: PipelineConfig,
    programs_paths: ProgramsPaths,
    output_dir: Path,
    parsed_subcommand: clapy.ParsedCommand,
) -> None:
    if len(list(output_dir.iterdir())) > 0:
        logger.warning("The output directory for the analysis pipeline is not empty")

    with TemporaryDirectory(
        dir=output_dir, delete=logger.level != logging.DEBUG
    ) as temp:
        temp_path = Path(temp)
        input_files: list[Path] = []
        assembly_path: Path

        match mode:
            case AnalysisMode.READS:
                input_files.extend(
                    list(parsed_subcommand.get_many("reads-paths", force=True))
                )
                assembly_results: Path | None = reads_assembly(
                    logger,
                    temp_path,
                    config,
                    programs_paths,
                    parsed_subcommand.get_one("query-path"),
                    output_dir,
                    list(parsed_subcommand.get_many("reads-paths", force=True)),
                    parsed_subcommand.get_flag("preprocessing"),
                    parsed_subcommand.get_flag("paired-ends"),
                )

                if assembly_results is None:
                    logger.critical("Couldn't assemble reads, aborting pipeline")
                    sys.exit(1)

                assembly_path = assembly_results

            case AnalysisMode.ASSEMBLY:
                # Check input assembly
                if not parsed_subcommand.get_one("assembly-path").is_file():
                    logger.critical(
                        f"Assembly file '{parsed_subcommand.get_one('assembly-path')}' not found, aborting pipeline."
                    )
                    sys.exit(1)
                else:
                    assembly_path = cast(
                        Path, parsed_subcommand.get_one("assembly-path")
                    )
                    input_files.append(assembly_path)

        # ====== Identification ======
        identification_results_path: Path = ensure_dir(
            logger,
            output_dir / config.identification.sub_dir,
            True,
            "can't run identification step. Aborting analysis pipeline",
        )

        identification_results: IdentificationResults
        with console.status("Identifying assembly features...", spinner="dots"):
            identification_results: IdentificationResults | None = run_identification(
                logger,
                config.identification,
                programs_paths,
                temp_path,
                identification_results_path,
                parsed_subcommand.get_one("query-path"),
                assembly_path,
            )

            if identification_results is None:
                logger.critical(
                    "Couldn't execute identification step. Aborting analysis pipeline"
                )
                sys.exit(1)

        logger.info("Identification completed")

        # ====== Differences Analysis ======
        difference_analysis_results_path: Path = ensure_dir(
            logger,
            output_dir / config.differences_analysis.sub_dir,
            True,
            "can't run differences analysis step. Aborting analysis pipeline",
        )

        differences_analysis_results: list[DifferencesAnalysisResults]
        hits_per_feature: HitsPerFeature = {}
        with console.status(
            "Analysis of query and assembly features differences...", spinner="dots"
        ):
            results: tuple[list[DifferencesAnalysisResults], HitsPerFeature] | None = (
                run_differences_analysis(
                    logger,
                    config.differences_analysis,
                    programs_paths,
                    temp_path,
                    difference_analysis_results_path,
                    identification_results["hits"],
                    assembly_path,
                    parsed_subcommand.get_one("query-path"),
                )
            )

            if results is None:
                logger.critical(
                    "No features found for differences analysis. Aborting analysis pipeline"
                )
                sys.exit(1)

            differences_analysis_results, hits_per_feature = results

            succesful_differences_analysis_features: list[str] = list(
                map(lambda x: x["feature_id"], differences_analysis_results)
            )

        logger.info("Differences analysis completed")

        # Check for additional phage matches with online-mode or metadata
        metadata_files: list[Path] = list(
            cast(tuple[Path, ...], parsed_subcommand.get_many("metadata", force=True))
        )
        phage_targets: dict[str, set[str]] = {}
        # Unimplemented online mode code
        # if parsed_subcommand.get_flag("online-mode"):
        #     if len(metadata_files) > 0:
        #         logger.info("Executing online mode analysis, ignoring metadata")
        if len(metadata_files) > 0:
            phage_targets: dict[str, set[str]] = metadata_analysis(
                logger,
                set(map(lambda x: x["feature_id"], identification_results["hits"])),
                metadata_files,
            )

        with console.status("Generating analysis report...", spinner="noise"):
            generate_report(
                AnalysisResults(
                    mode=mode,
                    input_files=input_files,
                    hits_per_feature=hits_per_feature,
                    phage_targets=phage_targets,
                ),
                config.report.theme,
                succesful_differences_analysis_features,
                config.report.nucleotide_database,
                output_dir,
            )

        logger.info(
            f"Analysis finished successfully! Results stored in '{output_dir}'",
        )
