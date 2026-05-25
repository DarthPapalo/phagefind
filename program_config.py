import logging
import sys
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import dacite

DEFAULT_PROGRAMS_PATHS_FILE: Path = Path(__file__).parent / "programs_paths.toml"
DEFAULT_CONFIG_FILE: Path = Path(__file__).parent / "config.toml"

logger: logging.Logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class ProgramsPaths:
    samtools: Path
    fastp: Path
    spades: Path

    @dataclass(frozen=True, slots=True)
    class BlastplusPaths:
        makeblastdb: Path
        blastn: Path

    blastplus: BlastplusPaths

    @dataclass(frozen=True, slots=True)
    class MummerPaths:
        nucmer: Path
        show_aligns: Path
        dnadiff: Path

    mummer: MummerPaths

    def load(programs_path: Path) -> ProgramsPaths:
        with programs_path.open("rb") as f:
            raw_config: dict[str, Any] = tomllib.load(f)
            return dacite.from_dict(
                ProgramsPaths, raw_config, dacite.Config(type_hooks={Path: Path})
            )


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    """
    Pipeline config class with default values.
    """

    @dataclass
    class PreprocessingConfig:
        sub_dir: str = "preprocessing"
        method: Literal["fastp"] = "fastp"

        @dataclass
        class PreprocessingFiltersConfig:
            n_base_limit: int = 5
            minimum_average_quality: int = 0
            qualified_phred_quality: int = 15
            unqualified_max_percentage: int = 40
            minimum_length: int = 75
            maximum_length: int = 350

        filters: PreprocessingFiltersConfig = field(
            default_factory=PreprocessingFiltersConfig
        )

    preprocessing: PreprocessingConfig = field(default_factory=PreprocessingConfig)

    @dataclass
    class AssemblyConfig:
        sub_dir: str = "assembly"
        method: Literal["SPAdes"] = "SPAdes"
        use_file: Literal["scaffolds", "contigs"] = "scaffolds"

    assembly: AssemblyConfig = field(default_factory=AssemblyConfig)

    @dataclass
    class IdentificationConfig:
        sub_dir: str = "identification"
        method: Literal["BLAST"] = "BLAST"

        @dataclass
        class IdentificationFiltersConfig:
            enable_dust: bool = False
            min_query_coverage_percentage: int = 5
            max_evalue: float = 1e-10

        filters: IdentificationFiltersConfig = field(
            default_factory=IdentificationFiltersConfig
        )

    identification: IdentificationConfig = field(default_factory=IdentificationConfig)

    @dataclass
    class DifferencesAnalysisConfig:
        sub_dir: str = "differences"
        generate_alignments: bool = True

        @dataclass
        class DifferencesAnalysisOptionsConfig:
            break_length: int = 200
            minimum_cluster: int = 65
            maximum_gap: int = 90
            minimum_match: int = 20

        options: DifferencesAnalysisOptionsConfig = field(
            default_factory=DifferencesAnalysisOptionsConfig
        )

    differences_analysis: DifferencesAnalysisConfig = field(
        default_factory=DifferencesAnalysisConfig
    )

    @dataclass
    class VisualizationConfig:
        standalone_graphs: bool = True
        color_palette: Literal["custom", "okabeito"] = "okabeito"

    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)

    @dataclass
    class OnlineModeConfig:
        email: str = ""

    online_mode: OnlineModeConfig = field(default_factory=OnlineModeConfig)

    @dataclass
    class ReportConfig:
        theme: Literal["light", "dark"] = "light"
        nucleotide_database: Literal["ENA", "NCBI"] = "NCBI"

    report: ReportConfig = field(default_factory=ReportConfig)

    def load(config_path: Path) -> PipelineConfig:
        """
        Loads a config TOML file path into a PipelineConfig dataclass.
        """
        with config_path.open("rb") as f:
            raw_config: dict[str, Any] = tomllib.load(f)
            return dacite.from_dict(PipelineConfig, raw_config)


def load_config(config_path: Path) -> PipelineConfig:
    return (
        PipelineConfig.load(config_path) if config_path.exists() else PipelineConfig()
    )


def load_programs_paths(programs_paths: Path) -> ProgramsPaths:
    if programs_paths.exists():
        return ProgramsPaths.load(programs_paths)
    else:
        logger.critical(f"Can't load programs paths file at '{programs_paths}'")
        sys.exit(1)
