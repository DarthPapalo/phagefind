import logging
import re
import ssl
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterator

import clapy
from rich.progress import BarColumn, MofNCompleteColumn, Progress, TextColumn, track

from program_config import ProgramsPaths

DOWNLOAD_DATA_COMMAND_NAME = "download-data"

AVAILABLE_DATA_SOURCES: tuple[str, ...] = (
    "RefSeq",
    "Genbank",
    "DDBJ",
    "EMBL",
    "PhagesDB",
    "GPD",
    "GVD",
    "MGV",
    "TemPhD",
    "CHVD",
    "IGVD",
    "IMG_VR",
    "GOV2",
    "STV",
)

AVAILABLE_DATA_TYPES: tuple[str, ...] = ("Genome", "Gene", "Metadata")

# fmt: off
cli: clapy.Command = (
    clapy.Command(DOWNLOAD_DATA_COMMAND_NAME)
        .help("Download bacteriophage data from the PhageScope online database.")
        .arguments([
                clapy.Arg("data-sources")
                    .short("-s")
                    .long("--sources")
                    .nargs(range(1, len(AVAILABLE_DATA_SOURCES)))
                    .help("Data sources that will be used to generate the database.")
                    .default(("RefSeq", "Genbank", "EMBL", "PhagesDB"))
                    .valid_values(AVAILABLE_DATA_SOURCES),
                clapy.Arg("data-types")
                    .short("-t")
                    .long("--data-types")
                    .nargs(range(1, len(AVAILABLE_DATA_TYPES)))
                    .help("Data types that can be selected to generate the database.")
                    .default(("Genome", "Metadata"))
                    .value_parser(str.capitalize)
                    .valid_values(AVAILABLE_DATA_TYPES),
                clapy.Arg("keep-individual")
                    .long("--keep-individual")
                    .help(
                        "Keep individual FASTA and GFF3 files per source. By default, these are deleted after merging."
                    )
                    .flag(),
                clapy.Arg("verify")
                    .long("--no-verify")
                    .flag(False)
                    .help("Don't check for SSL certificates before downloading data."),
            ])
)
# fmt: on


def fasta_from_gff(
    logger: logging.Logger,
    samtools_path: Path,
    fasta_path: Path,
    gff_path: Path,
    target_feature: str,
    result_path: Path,
) -> None:
    """
    Extracts the sequences from a given feature from a Fasta file to another fasta file using a gff3 annotation file.
    """
    with TemporaryDirectory(
        dir=result_path.parent, delete=logger.level != logging.DEBUG
    ) as work_dir:
        work_dir_path: Path = Path(work_dir)
        regions_temp: Path = work_dir_path / "regions.tmp"
        headers_temp: Path = work_dir_path / "headers.tmp"

        with (
            regions_temp.open("w") as regions_file,
            headers_temp.open("w") as headers_file,
        ):
            for line in gff_path.read_text().splitlines():
                columns: list[str] = line.split("\t")
                if len(columns) < 9:
                    continue
                if columns[2] == target_feature:
                    match: re.Match[str] | None = re.search(r"ID=([^;]*);", columns[8])
                    if match is None:
                        logger.error(
                            "No sequence ID found for targeted feature gff annotation, skipping..."
                        )
                        continue
                    id: str = match.group(1)

                    regions_file.write(f"{columns[0]}:{columns[3]}-{columns[4]}\n")
                    headers_file.write(f">{columns[0]}|{target_feature}:{id}\n")

        # Duplicate fasta sequences as some genes are anotated for circular DNA
        duplicated_fasta_temp: Path = work_dir_path / "duplicated_fasta.tmp"
        with duplicated_fasta_temp.open("w") as duplicated_fasta_file:
            for line in fasta_path.read_text().splitlines():
                if line.startswith(">"):
                    duplicated_fasta_file.write(line + "\n")
                else:
                    duplicated_fasta_file.write(line * 2 + "\n")

        regions_temp: Path = work_dir_path / "regions.tmp"
        samtools_regions_args: list[str] = [
            str(samtools_path),
            "faidx",
            str(duplicated_fasta_temp),
            "--region-file",
            str(regions_temp),
        ]

        with regions_temp.open("w") as regions_file:
            try:
                subprocess.run(samtools_regions_args, stdout=regions_file, check=True)
            except subprocess.CalledProcessError:
                logger.error(
                    "Couldn't extract regions with samtools in 'fasta_from_gff'"
                )
                return

        headers: Iterator[str] = iter(headers_temp.read_text().splitlines())

        with regions_temp.open() as regions_file, open(result_path, "w") as result_file:
            for line in regions_file:
                if line.startswith(">"):
                    result_file.write(next(headers) + "\n")
                else:
                    result_file.write(line)


def get_data(
    logger: logging.Logger,
    loaded_programs_paths: ProgramsPaths,
    output_path: Path,
    data_sources: list[str],
    data_types: list[str],
    keep_individual: bool,
    verify_ssl: bool,
):
    """
    Generate databases from online resources that can be used for the pipeline.
    """
    lower_data_types: set[str] = set(map(str.lower, data_types))

    ssl_context: ssl.SSLContext = ssl.create_default_context()
    if not verify_ssl:
        ssl_context.check_hostname = False
        ssl_context.verify_mode = ssl.VerifyMode.CERT_NONE
    try:
        if "metadata" in lower_data_types:
            meta_dir: Path = output_path / "meta"
            meta_dir.mkdir(parents=True)

            with Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
            ) as p:
                for source in p.track(
                    data_sources,
                    description="Downloading metadata from sources",
                ):
                    # ====== Download Metadata files ======
                    meta_url = f"https://phageapi.deepomics.org/files/Download/Phage_meta_data/{source.lower()}_phage_meta_data.tsv"
                    meta_output: Path = meta_dir / f"{source}_phage_metadata.tsv"
                    with urllib.request.urlopen(
                        meta_url, context=ssl_context
                    ) as response:
                        with open(meta_output, "wb") as f:
                            while True:
                                chunk = response.read(8192)
                                if not chunk:
                                    break
                                f.write(chunk)

            logger.info("Successfully acquired Metadata")

        if "genome" in lower_data_types:
            sequences_dir: Path = output_path / "sequences"
            sequences_dir.mkdir(parents=True)

            all_genomes_path: Path = sequences_dir / "all_genomes.fasta"

            for source in track(
                data_sources,
                description="Downloading genome data from sources",
                show_speed=False,
            ):
                # ====== Download tar.gz Sequence files ======
                tar_url = f"https://phageapi.deepomics.org/download/phage/fasta/?datasource={source}"
                temp_tar: Path = sequences_dir / f"tmp-{source}.tar.gz"
                with urllib.request.urlopen(tar_url, context=ssl_context) as response:
                    with open(temp_tar, "wb") as f:
                        while True:
                            chunk = response.read(8192)
                            if not chunk:
                                break
                            f.write(chunk)

                # ====== Extract ======
                source_dir: Path = sequences_dir / source
                source_dir.mkdir(parents=True, exist_ok=True)

                with tarfile.open(temp_tar, "r:gz") as tar:
                    tar.extractall(path=source_dir)

                temp_tar.unlink()

                # ====== Join FASTA files ======
                sequences_output_path: Path = sequences_dir / f"{source}_genomes.fasta"
                with sequences_output_path.open("ab") as sequences_output_file:
                    for fasta_path in source_dir.rglob("*.fasta"):
                        with fasta_path.open("rb") as fasta_file:
                            sequences_output_file.write(fasta_file.read())

                # ====== Erase Individual files ======
                if not keep_individual:
                    for p in source_dir.rglob("*.fasta"):
                        if p.is_file():
                            p.unlink()
                    for p in source_dir.rglob("*"):
                        if p.is_dir():
                            p.rmdir()
                    source_dir.rmdir()

            # Join all data sources sequences into one
            genomes_fasta_paths: list[Path] = list(
                sequences_dir.glob("*_genomes.fasta")
            )
            with all_genomes_path.open("w") as all_genomes_file:
                for genome_path in genomes_fasta_paths:
                    all_genomes_file.write(genome_path.read_text() + "\n")
                    if not keep_individual:
                        genome_path.unlink()

            logger.info("Successfully acquired Genome data")

            if "gene" in lower_data_types:
                annotations_dir = output_path / "annotation"
                annotations_dir.mkdir()

                for source in track(
                    data_sources,
                    description="Downloading annotation data from sources",
                    show_speed=False,
                ):
                    annotation_url = f"https://phageapi.deepomics.org/fasta/phage_sequence/phage_gff3/{source}.gff3"
                    annotation_path: Path = annotations_dir / f"{source}.gff3"
                    with urllib.request.urlopen(
                        annotation_url, context=ssl_context
                    ) as response:
                        with open(annotation_path, "wb") as f:
                            while True:
                                chunk = response.read(8192)
                                if not chunk:
                                    break
                                f.write(chunk)

                merged_annotations_path: Path = annotations_dir / "all_annotations.gff3"

                annotations_paths: list[Path] = list(annotations_dir.glob("*.gff3"))
                with merged_annotations_path.open("w") as merged_annotations_file:
                    for annotation_path in annotations_paths:
                        merged_annotations_file.write(
                            annotation_path.read_text() + "\n"
                        )
                        if not keep_individual:
                            annotation_path.unlink()

                fasta_from_gff(
                    logger,
                    loaded_programs_paths.samtools,
                    all_genomes_path,
                    merged_annotations_path,
                    "gene",
                    sequences_dir / "all_genes.fasta",
                )

                logger.info("Successfully acquired Gene data")

        elif "gene" in lower_data_types:
            logger.error(
                "Error downloading data: Genome data is required to extract gene data"
            )
            sys.exit(1)

    except urllib.error.URLError as ue:
        logger.critical(f"Error downloading data: {ue.reason}")
        sys.exit(1)
    except FileExistsError as fe:
        logger.critical(
            f"Error downloading data: File or directory already exists: '{fe.filename}'"
        )
        sys.exit(1)
