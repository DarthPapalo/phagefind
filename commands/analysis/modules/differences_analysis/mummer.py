import logging
import subprocess
from pathlib import Path
from typing import TypedDict

from program_config import PipelineConfig


class MummerResults(TypedDict):
    alignments_path: Path
    snps_path: Path
    report_path: Path
    one_to_one_coords_path: Path


def run_mummer(
    logger: logging.Logger,
    nucmer_path: Path,
    show_aligns_path: Path,
    dnadiff_path: Path,
    config: PipelineConfig.DifferencesAnalysisConfig.DifferencesAnalysisOptionsConfig,
    reference_fasta_path: Path,
    reconstructed_fasta_path: Path,
    feature_id: str,
    reconstructed_id: str,
    workdir: Path,
) -> MummerResults | None:
    """
    Executes MUMMER for the differences analysis, returns None in case errors occurred.
    """
    # Mummer delta file generation
    delta_file: Path = workdir / "differences_analysis_mummer.delta"
    mmumer_args: list[str] = [
        str(nucmer_path),
        "--delta",
        str(delta_file),
        "--breaklen",
        str(
            config.break_length
        ),  # Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (200)
        "--mincluster",
        str(
            config.minimum_cluster
        ),  # Sets the minimum length of a cluster of matches (65)
        "--maxgap",
        str(
            config.maximum_gap
        ),  # Set the maximum gap between two adjacent matches in a cluster (90)
        "--minmatch",
        str(
            config.minimum_match
        ),  # Set the minimum length of a single exact match (20)
        str(reference_fasta_path),
        str(
            reconstructed_fasta_path
        ),  # Reference - Query // In our case we are using the reconstructed sequence as the query !!
    ]
    try:
        subprocess.run(
            mmumer_args, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL, check=True
        )
    except subprocess.CalledProcessError:
        logger.warning(
            f"Mummer couldn't generate the delta file for {reference_fasta_path} / {reconstructed_fasta_path}"
        )
        return None

    # Show alignments
    alignments_file: Path = workdir / "alignments_mummer.ali"
    showalign_args: list[str] = [
        str(show_aligns_path),
        "-r",
        str(delta_file),
        feature_id,
        reconstructed_id,
    ]

    with alignments_file.open("w") as alignments_output:
        try:
            subprocess.run(
                showalign_args,
                stderr=subprocess.DEVNULL,
                stdout=alignments_output,
                check=True,
            )
        except subprocess.CalledProcessError:
            logger.warning(
                f"MUMMER couldn't generate alignments for feature '{feature_id}'"
            )
            return None

    # Dna differenes, generates snps file
    dna_diff_dir: Path = workdir / "dna_diff"
    dna_diff_dir.mkdir()

    dnadiff_args: list[str] = [
        str(dnadiff_path),
        "-d",
        str(delta_file),
        "-p",
        str(dna_diff_dir / "out"),
    ]

    # This command prints redundant information to stderr...
    try:
        subprocess.run(
            dnadiff_args,
            stderr=subprocess.STDOUT,
            stdout=subprocess.DEVNULL,
            check=True,
        )
    except subprocess.CalledProcessError:
        logger.warning(
            f"Mummer couldn't run the DNA difference analysis for '{delta_file}' delta file"
        )
        return None

    # Clean SNPs file
    clean_snps: Path = workdir / "snps.tsv"
    with (
        (dna_diff_dir / "out.snps").open() as mummer_snps_file,
        clean_snps.open("w") as new_snps_file,
    ):
        new_snps_file.write(
            "SNP position in reference phage seq\t"
            "Nt in reference phage seq\t"
            "Nt in found phage seq\t"
            "SNP position in found phage seq\t"
            "Distance from this SNP to the nearest mismatch\t"
            "Distance from this SNP to nearest sequence end\t"
            "Length of reference phage seq\t"
            "Length of found phage seq\t"
            "Reference phage seq strand\t"
            "Found phage seq strand\t"
            "Phage reference SeqID\t"
            "Found phage SeqID\t"
            "\n"
        )
        new_snps_file.write(mummer_snps_file.read())

    return MummerResults(
        alignments_path=alignments_file,
        snps_path=clean_snps,
        report_path=dna_diff_dir / "out.report",
        one_to_one_coords_path=dna_diff_dir / "out.1coords",
    )
