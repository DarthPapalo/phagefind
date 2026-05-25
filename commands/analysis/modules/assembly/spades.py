import logging
import subprocess
from pathlib import Path

from program_config import PipelineConfig


def run_spades(
    logger: logging.Logger,
    spades_path: Path,
    config: PipelineConfig.AssemblyConfig,
    reads_paths: list[Path],
    paired_ends: bool,
    workdir: Path,
) -> tuple[Path | None, Path | None, Path]:
    """
    Returns a tuple with the (contigs, scaffolds, log) files.
    """

    def check_all_fastq(reads: list[Path]) -> bool:
        """
        Checks if all reads are in fastq format.
        """
        for r in reads:
            if "".join(r.suffixes) not in [".fq", ".fastq", ".fq.gz", ".fastq.gz"]:
                return False
        return True

    spades_args: list[str] = [
        str(spades_path),
        "--careful",
        "-o",
        str(workdir),
    ]

    if not check_all_fastq(reads_paths):
        # If not all files are fastq we can only perform the assembly, no reads error correction
        spades_args.append("--only-assembler")

    if paired_ends:
        for i, file in enumerate(reads_paths):
            assert file.exists()
            pair_idx = i // 2 + 1
            read_idx = i % 2 + 1
            spades_args.extend([f"--pe{pair_idx}-{read_idx}", str(file)])
    else:
        for file in reads_paths:
            assert file.exists()
            spades_args.extend(["-s", str(file)])

    log_path: Path = workdir / "SPAdes.log"

    try:
        subprocess.run(
            spades_args, stderr=subprocess.STDOUT, stdout=log_path.open("a"), check=True
        )
    except subprocess.CalledProcessError:
        logger.error("SPAdes failed to assemble a genome from the reads")
        return (None, None, log_path)

    contigs_path = workdir / "contigs.fasta"
    scaffolds_path = workdir / "scaffolds.fasta"

    return (contigs_path, scaffolds_path, log_path)
