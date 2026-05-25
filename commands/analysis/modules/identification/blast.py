import logging
import subprocess
from pathlib import Path

from program_config import PipelineConfig


def run_blast(
    logger: logging.Logger,
    makeblastdb_path: Path,
    blastn_path: Path,
    config: PipelineConfig.IdentificationConfig.IdentificationFiltersConfig,
    query_path: Path,
    genome_path: Path,
    work_dir: Path,
) -> Path | None:
    """
    Runs the makeblastdb (To generate a blast database) and the blastn commands
    Return None in case of errors
    """
    blastdb_path: str = str(work_dir / "assembly_index") + "/db"

    # Index assembly for blast
    # May truncate sequence names with spaces TODO
    makeblastsb_args: list[str] = [
        str(makeblastdb_path),
        "-in",
        str(genome_path),
        "-input_type",
        "fasta",
        "-parse_seqids",
        "-dbtype",
        "nucl",
        "-out",
        str(blastdb_path),
    ]

    try:
        subprocess.run(
            makeblastsb_args,
            stderr=subprocess.STDOUT,
            stdout=subprocess.DEVNULL,
            check=True,
        )
    except subprocess.CalledProcessError:
        logger.error(f"BLAST+ failed to make the db for '{genome_path}'")
        return None

    blast_result_path: Path = work_dir / "blast_hits.tsv"

    # Executing BLASTn
    blast_args: list[str] = [
        str(blastn_path),
        "-query",
        str(query_path),
        "-db",
        str(blastdb_path),
        "-dust",
        "yes" if config.enable_dust else "no",
        "-evalue",
        str(config.max_evalue),
        "-qcov_hsp_perc",
        str(config.min_query_coverage_percentage),
        "-outfmt",
        "6 qseqid sseqid qstart qend sstart send frames qcovhsp qcovs evalue bitscore score sseq",
        "-out",
        str(blast_result_path),
    ]
    try:
        subprocess.run(
            blast_args, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL, check=True
        )
    except subprocess.CalledProcessError:
        logger.error(f"BLAST+ failed to execute blastn for query '{query_path}'")
        return None

    return blast_result_path
