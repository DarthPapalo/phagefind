import logging
from pathlib import Path

from program_config import PipelineConfig, ProgramsPaths

from .fastp import run_fastp


def run(
    logger: logging.Logger,
    config: PipelineConfig.PreprocessingConfig,
    programs_paths: ProgramsPaths,
    output_path: Path,
    reads_paths: list[Path],
    paired_ends: bool,
) -> list[Path] | None:

    fastq_reads: list[Path] = []
    fasta_reads: list[Path] = []
    for r in reads_paths:
        if "".join(r.suffixes) in [".fq", ".fastq", ".fq.gz", ".fastq.gz"]:
            fastq_reads.append(r)
        else:
            fasta_reads.append(r)

    if len(fastq_reads) == 0:
        logger.warning(
            "No reads in fastq format found, skipping preprocessing of reads"
        )
        return reads_paths

    if config.method == "fastp":
        preprocessed_reads = run_fastp(
            logger,
            programs_paths.fastp,
            config.filters,
            fastq_reads,
            paired_ends,
            output_path,
        )
        if preprocessed_reads is None:
            return None

        # Check preprocessed reads size
        for f in preprocessed_reads:
            if f.stat().st_size == 0:
                logger.error(
                    f"Preprocessing of reads ended with no valid reads for '{f}'. Check your preprocessing configuration parameters"
                )
                return None

        return fasta_reads + preprocessed_reads
    else:
        logger.error(f"Unknown preprocessing method: {config.method}")
        return None
