import logging
import sys
from pathlib import Path

import clapy
from rich.progress import track

from _commons import console
from commands.analysis.modules.differences_analysis import (
    RESULTS_ONE_TO_ONE_COORDS_FILE,
    RESULTS_SNPS_FILE,
)
from program_config import PipelineConfig

from .graph_generation import AVAILABLE_GRAPHS, generate_graphs

VISUALIZATION_COMMAND_NAME = "visualization"

ALL_FEATURES_TOKEN: str = "ALL"

# fmt: off
cli: clapy.Command = (
    clapy.Command(VISUALIZATION_COMMAND_NAME)
        .help("Generate different visualisations of the analysis results.")
        .arguments([
            clapy.Arg("features-dir")
                .help("Directory with the differences analysis results for the different features identified.")
                .long("--features-dir")
                .value_parser(Path),
            clapy.Arg("features-ids")
                .help(f"The ID(s) of the feature(s) to generate visualziations for (Use '{ALL_FEATURES_TOKEN}' for all the features).")
                .nargs("+"),
            clapy.Arg("graphs")
                .help("The graphs you want to generate (By default all are generated).")
                .long("--graphs")
                .nargs(range(1, len(AVAILABLE_GRAPHS)))
                .valid_values(AVAILABLE_GRAPHS)
                .default(AVAILABLE_GRAPHS),
        ])
)
# fmt: on


def run(
    logger: logging.Logger,
    config: PipelineConfig.VisualizationConfig,
    export_dir: Path,
    features_dir: Path,
    features_ids: tuple[str, ...],
    graphs: tuple[str, ...],
) -> None:
    def path_last_part(p: Path) -> str:
        return str(p.parts[-1])

    features_dirs: list[Path] = []
    try:
        features_dirs: list[Path] = list(features_dir.iterdir())
    except FileNotFoundError:
        logger.critical(
            "Features directory is empty, aborting visualizations generation"
        )
        sys.exit(1)

    if ALL_FEATURES_TOKEN in features_ids:
        # Select all features (directories) if '*' is used as a feature ID.
        features_ids = tuple(map(path_last_part, filter(Path.is_dir, features_dirs)))

    for feature in track(features_ids, "Generating visualizations..."):
        if feature not in list(map(path_last_part, features_dirs)):
            logger.warning(
                f"Feature '{feature}' not found in features dir, skipping..."
            )
            continue

        feature_export_dir: Path = export_dir / feature
        if feature_export_dir.exists():
            logger.error(
                f"Visualizations output path already has a dir for feature '{feature}', skipping..."
            )
            continue

        with console.status(
            f"Generating visualizations for feature [cyan]{feature}[/cyan]"
        ):
            generate_graphs(
                logger,
                config,
                graphs,
                (features_dir / feature / RESULTS_SNPS_FILE),
                (features_dir / feature / RESULTS_ONE_TO_ONE_COORDS_FILE),
                feature_export_dir,
            )

    logger.info(
        f"Visualizations generated successfully! Results stored in '{export_dir}'"
    )
