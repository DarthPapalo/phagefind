import logging
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import TypedDict

from commands.analysis.modules.identification import Hit
from program_config import PipelineConfig, ProgramsPaths

from .mummer import MummerResults, run_mummer

RESULTS_REPORT_FILE = "differences_summary_report.txt"
RESULTS_SNPS_FILE = "snps.tsv"
RESULTS_ALIGNMENT_FILE = "alignments.txt"
RESULTS_ONE_TO_ONE_COORDS_FILE = "one_to_one_coords.1coords"
RESULTS_RECONSTRUCTED_FILE = "reconstructed.fasta"

type HitsPerFeature = dict[str, list[Hit]]


class DifferencesAnalysisResults(TypedDict):
    feature_id: str
    report_path: Path
    snps_path: Path
    alignments_path: Path
    one_to_one_coords_path: Path
    reconstructed_path: Path
    hits: list[Hit]


def extract_sequence(
    samtools_path: Path,
    output_path: Path,
    fasta_path: Path,
    seqid: str,
    start_end: tuple[int, int] | None = None,
    reverse_complement: bool = False,
) -> Path | None:
    """
    Appends a sequence from a Fasta file using its seqid into its own file.
    Returns the new file path.
    Uses samtools.
    """
    extract_seq_args: list[str] = [
        str(samtools_path),
        "faidx",
        str(fasta_path),
        f"{seqid}{f':{start_end[0]}-{start_end[1]}' if start_end is not None else ''}",
        "--length",
        "0",
    ]
    if reverse_complement:
        extract_seq_args.extend(["--reverse-complement", "--mark-strand", "no"])

    try:
        with open(output_path, "a") as output_file:
            subprocess.run(extract_seq_args, stdout=output_file, check=True)
    except subprocess.CalledProcessError:
        return None

    return output_path


def parse_location(loc_str: str) -> tuple[int, int]:
    """
    Parses a location string such as '<start>-<end>' into a tuple[int, int].
    """
    start, end = tuple(map(int, loc_str.split("-")))
    return (min(start, end), max(start, end))


def get_interest_hits(
    hits: list[Hit], min_cover: int = 150
) -> list[tuple[int, tuple[int, int]]]:
    """
    Return the hits indexes required for the maximum feature cover possible together
    with the feature region for each one. list[(hit-idx, (feature-start, feature-end))]
    """
    # If equal feature_loc in different hits use the one with the best Bitscore
    # scores: dict[(start, end), (index, bitscore)]
    scores: dict[tuple[int, int], tuple[int, int]] = defaultdict(lambda: (0, 0))
    for i, (start, end), bitscore in [
        (i, parse_location(hit["feature_loc"]), hit["bitscore"])
        for i, hit in enumerate(hits)
    ]:
        if bitscore > scores[(start, end)][1]:
            scores[(start, end)] = (i, bitscore)

    # intervals: list[(hit_idx, (feature_start, feature_end))]
    intervals: list[tuple[int, tuple[int, int]]] = [
        (scores[(start, end)][0], (start, end)) for (start, end) in scores.keys()
    ]
    intervals.sort(key=lambda x: x[1][0])  # Sorted by interval start
    selected: list[tuple[int, tuple[int, int]]] = [intervals[0]]

    # Select larger starter interval
    for i in intervals[1:]:
        if i[1][1] > selected[-1][1][1]:  # Extends what is already covered:
            if (
                i[1][0] > selected[-1][1][1]
            ):  # Starts after what is already covered -> Empty space generated
                selected.append((i[0], i[1]))
            elif (
                selected[-1][1][1] + 1 < i[1][1] - min_cover
            ):  # Starts on already covered location AND at leasts covers min_cover
                selected.append((i[0], (selected[-1][1][1] + 1, i[1][1])))

    return selected


def reconstruct_sequence(
    hits: list[Hit],
    output_path: Path,
) -> str:
    """
    Writes the reconstruction of found feature using a list of BLAST hits
    that have the highest coverage possible combined.
    Returns the name of the SeqID used in the created FASTA file.
    """
    # We get the hits that return the maximum reference coverage
    # This hits contain the bacterial genome regions of the parts that are similar to the feature
    interest_hits: list[tuple[int, tuple[int, int]]] = get_interest_hits(hits)

    with output_path.open("a") as f:
        last_pos = 0
        for i in interest_hits:
            current_hit: Hit = hits[i[0]]
            # Add possible gaps in the genome between hits
            previous_gap = "-" * (i[1][0] - last_pos)
            # Dont add gaps before first hit
            if last_pos != 0 and len(previous_gap) > 0:
                f.write(previous_gap)
            f.write(current_hit["aligned_sequence"])
            last_pos = i[1][1]

    # Join all the BLAST aligned bacteria genome fragments to create the "similar-to-feature" sequence
    full_sequence: str = ""
    with open(output_path) as dirty:
        full_sequence = dirty.read().replace("\n", "")

    reconstructed_seqid = f"ReconstructedFeature|similar_to={hits[0]['feature_id']}"
    with open(output_path, "w") as clean:
        clean.write(">" + reconstructed_seqid + "\n")
        # 70 as a limit of characters for the FASTA file
        BASES_PER_LINE = 70
        lines: list[str] = [
            full_sequence[i : i + BASES_PER_LINE]
            for i in range(0, len(full_sequence), BASES_PER_LINE)
        ]
        clean.write("\n".join(lines))

    return reconstructed_seqid


def run(
    logger: logging.Logger,
    config: PipelineConfig.DifferencesAnalysisConfig,
    programs_paths: ProgramsPaths,
    temp_path: Path,
    results_path: Path,
    hits_list: list[Hit],
    assembly_path: Path,
    query_path: Path,
) -> tuple[list[DifferencesAnalysisResults], HitsPerFeature] | None:
    """
    Runs the a differences analysis for the features detected in the bacteria assembly
    compared to the reference sequences. The last step of the pipeline.
    """

    results: list[DifferencesAnalysisResults] = []
    hits_per_feature: dict[str, list[Hit]] = defaultdict(list)

    # We group the hits from the identification step
    for hit in hits_list:
        hits_per_feature[hit["feature_id"]].append(hit)

    for feature in hits_per_feature.keys():

        def sanitize_dirname(name: str, replacement: str = "_") -> str:
            """Replace characters invalid in directory names with a safe alternative."""
            import re

            INVALID_CHARS: str = r"[<>:\"/\\|?*\x00-\x1f]"
            sanitized: str = re.sub(INVALID_CHARS, replacement, name)
            sanitized = sanitized.strip(". ")
            if not sanitized:
                sanitized = replacement
            return sanitized

        feature_results_path: Path = results_path / sanitize_dirname(feature)

        if (
            feature_results_path.is_dir()
            and len(list(feature_results_path.iterdir())) > 0
        ):
            logger.error(
                f"Can't run differences analysis for feature '{feature}', subdirectory exists and is not empty: '{results_path}'. Skipping feature..."
            )
            continue

        feature_workdir: Path = temp_path / "differences_analysis" / feature
        feature_workdir.mkdir(parents=True)

        # Extract individual query sequence incase its in a Multi-FASTA file
        query_fasta_path: Path | None = extract_sequence(
            programs_paths.samtools,
            feature_workdir / "query.fasta",
            query_path,
            feature,
        )

        if query_fasta_path is None:
            logger.error(
                f"Failed to extract query sequence from {query_path} with samtools, skipping feature"
            )
            continue

        # Create a feature reconstruction with the aligned BLAST hits from the bacteria genome
        reconstructed_fasta_path: Path = feature_workdir / "reconstructed.fasta"
        reconstructed_seqid: str = reconstruct_sequence(
            hits_per_feature[feature],
            reconstructed_fasta_path,
        )

        # Run mummer analysis for the feature
        analysis_results: MummerResults | None = run_mummer(
            logger,
            programs_paths.mummer.nucmer,
            programs_paths.mummer.show_aligns,
            programs_paths.mummer.dnadiff,
            config.options,
            query_fasta_path,
            reconstructed_fasta_path,
            feature,
            reconstructed_seqid,
            feature_workdir,
        )

        if analysis_results is None:
            logger.warning(
                f"Failed to execute the differences analysis for feature '{feature}', skipping..."
            )
            continue

        # Move interest files out of the temporary directory and store in results
        feature_results_path.mkdir(exist_ok=True)
        results.append(
            DifferencesAnalysisResults(
                feature_id=feature,
                report_path=analysis_results["report_path"].move(
                    feature_results_path / RESULTS_REPORT_FILE
                ),
                snps_path=analysis_results["snps_path"].move(
                    feature_results_path / RESULTS_SNPS_FILE
                ),
                alignments_path=analysis_results["alignments_path"].move(
                    feature_results_path / RESULTS_ALIGNMENT_FILE
                ),
                reconstructed_path=reconstructed_fasta_path.move(
                    feature_results_path / RESULTS_RECONSTRUCTED_FILE
                ),
                one_to_one_coords_path=analysis_results["one_to_one_coords_path"].move(
                    feature_results_path / RESULTS_ONE_TO_ONE_COORDS_FILE
                ),
                hits=hits_per_feature[feature],
            )
        )

    if len(results) > 0:
        return (results, hits_per_feature)
    else:
        return None
