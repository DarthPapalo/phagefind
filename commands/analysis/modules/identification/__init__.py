import logging
from pathlib import Path
from typing import TypedDict

from program_config import PipelineConfig, ProgramsPaths

from .blast import run_blast

ANNOTATION_FILE_NAME = "identification_annotation.tsv"


class Hit(TypedDict):
    feature_id: str
    genome_region: str
    feature_loc: str
    genome_region_loc: str
    bitscore: int
    evalue: float
    aligned_sequence: str


class IdentificationResults(TypedDict):
    """
    Results from the identification step of the pipeline.
    """

    annotation_path: Path
    hits: list[Hit]


def run(
    logger: logging.Logger,
    config: PipelineConfig.IdentificationConfig,
    programs_paths: ProgramsPaths,
    temp_path: Path,
    output_path: Path,
    query_path: Path,
    assembly_path: Path,
) -> IdentificationResults | None:
    """
    Returns an annotation TSV with the features found in the genome according to the query parameters.
    Used in the diferences analysis.
    """
    annotation_path: Path = output_path / ANNOTATION_FILE_NAME

    with open(annotation_path, "w") as clean_file:
        # Write header
        clean_file.write(
            "Identified feature\t"
            "Genome region\t"
            "Feature start-end\t"
            "Region start-end\t"
            "Query/Subject frames\t"
            "HSP coverage/Total feature coverage\t"
            "E-value\t"
            "Bitscore\t"
            "Raw score\t\n"
        )

        match config.method:
            case "BLAST":
                hits_file_path: Path | None = run_blast(
                    logger,
                    programs_paths.blastplus.makeblastdb,
                    programs_paths.blastplus.blastn,
                    config.filters,
                    query_path,
                    assembly_path,
                    temp_path,
                )

                if hits_file_path is None:
                    return None

                hits: list[Hit] = []
                # Write hits file with correct format
                with open(hits_file_path) as hits_file:
                    for line in list(
                        map(
                            lambda line: line.strip().split("\t"), hits_file.readlines()
                        )
                    ):
                        clean_file.write(
                            f"{line[0]}\t{line[1]}\t{line[2]}-{line[3]}\t{line[4]}-{line[5]}\t{line[6]}\t{line[7]}/{line[8]}\t{line[9]}\t{line[10]}\t{line[11]}\n"
                        )
                        hits.append(
                            Hit(
                                feature_id=line[0],
                                genome_region=line[1],
                                feature_loc=f"{line[2]}-{line[3]}",
                                genome_region_loc=f"{line[4]}-{line[5]}",
                                bitscore=int(line[7]),
                                evalue=float(line[9]),
                                aligned_sequence=line[12],
                            )
                        )
            case _:
                logger.error(f"Unknown identification method: {config.method}")
                return None

    return IdentificationResults(annotation_path=annotation_path, hits=hits)
