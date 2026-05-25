from enum import Enum
from pathlib import Path
from typing import TypedDict

from commands.analysis.modules.identification import Hit

METADATA_PHAGE_ACC_COLUMN = "Phage_ID"
METADATA_PHAGE_HOST_COLUMN = "Host"


class AnalysisMode(Enum):
    READS = 0
    ASSEMBLY = 1


class AnalysisResults(TypedDict):
    mode: AnalysisMode
    input_files: list[Path]
    hits_per_feature: dict[str, list[Hit]]
    phage_targets: dict[str, set[str]]
