import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Literal, cast

from pandas import DataFrame, read_table
from plotly.graph_objects import Figure

from program_config import PipelineConfig

AVAILABLE_GRAPHS: tuple[str, ...] = ("SNPs",)

UNALIGNED_GRAPH_VALUE = 1
MATCH_GRAPH_VALUE = 0
GAP_GRAPH_VALUE = -1
MISMATCH_GRAPH_VALUE = -2

# Color palettes
CUSTOM_PALETTE: list[str] = ["#009E73", "#D55E00", "#774cb0"]
OKABEITO_PALETTE: list[str] = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000",
]


@dataclass(slots=True)
class VisualizationResults:
    snps_graph_path: Path | None = None


def alignment_position_to_score(nts: tuple[str, str]) -> int:
    if nts[0] == nts[1]:
        return MATCH_GRAPH_VALUE
    elif any(x in nts for x in ("-", ".")):
        return GAP_GRAPH_VALUE
    else:
        return MISMATCH_GRAPH_VALUE


def _get_complete_segments(segments: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """
    From a list of segments returns a new list with contiguous segments merged and missing ones added.
    """
    segments: list[tuple[int, int]] = sorted(
        segments, key=lambda x: x[0]
    )  # Sorted by starting pos

    complete_segments: list[tuple[int, int]] = [segments[0]]

    for segment in segments[1:]:
        if (segment[0] - 1) == complete_segments[-1][
            1
        ]:  # Starts when last ends -> Merge
            merged_segment: tuple[int, int] = (complete_segments[-1][0], segment[1])
            complete_segments.pop(-1)
            complete_segments.append(merged_segment)
        else:  # Starts after the last one -> Add missing
            complete_segments.append((complete_segments[-1][1], segment[0]))  # Missing
            complete_segments.append(segment)

    return complete_segments


def generate_snps_graph(
    logger: logging.Logger,
    color_palette: list[str],
    snps_path: Path,
    one_to_one_coords_path: Path,
) -> Figure | None:

    snps_data: DataFrame = read_table(snps_path, sep="\t", header=0)
    if len(snps_data.index) <= 0:
        logger.error("Empty SNPs data")
        return None

    coords_data: DataFrame = read_table(one_to_one_coords_path, sep="\t", header=None)
    if len(snps_data.index) <= 0:
        logger.error("Empty one to one coords data")
        return None

    # Columns 2 and 3 (starting at index 0) represent start and end of the aligned *query* sequence
    query_regions_strands: list[Literal["+", "-"]] = [
        "+" if start < end else "-"
        for start, end in zip(coords_data.iloc[:, 2], coords_data.iloc[:, 3])
    ]
    aligned_query_regions: list[tuple[int, int]] = [
        (min(start, end), max(start, end))
        for start, end in zip(coords_data.iloc[:, 2], coords_data.iloc[:, 3])
    ]

    # Alignment start with final position (Contains gaps)
    total_alignment_range: tuple[int, int] = (
        min(pos for region in aligned_query_regions for pos in region),
        max(pos for region in aligned_query_regions for pos in region),
    )

    all_regions: list[tuple[int, int]] = _get_complete_segments(aligned_query_regions)

    # Positions in the found/reconstructed sequence
    # Pos: Score
    parsed_snps_scores: dict[int, int] = {
        pos: alignment_position_to_score((nt1, nt2))
        for pos, nt1, nt2 in zip(
            snps_data.iloc[:, 3], snps_data.iloc[:, 1], snps_data.iloc[:, 2]
        )
    }

    fig = Figure()

    for i, region in enumerate(all_regions):
        is_aligned_region: bool = i % 2 == 0
        line_df = DataFrame(
            dict(
                x=list(range(*region)),
                y=[
                    parsed_snps_scores[pos]
                    if pos in parsed_snps_scores
                    else MATCH_GRAPH_VALUE
                    for pos in range(*region)
                ]
                if is_aligned_region
                else [UNALIGNED_GRAPH_VALUE] * len(range(*region)),
            )
        )

        fig.add_scatter(
            x=line_df["x"],
            y=line_df["y"],
            mode="lines",
            line=dict(
                color=color_palette[0] if is_aligned_region else color_palette[1],
            ),
            showlegend=False,
            name="",
        )

    for i, (region, strand) in cast(
        Iterator[tuple[int, tuple[tuple[int, int], Literal["+", "-"]]]],
        enumerate(zip(aligned_query_regions, query_regions_strands)),
    ):
        if strand == "+":
            continue
        fig.add_vrect(
            x0=region[0],
            x1=region[1],
            fillcolor=color_palette[2],
            opacity=0.2,
            layer="below",
            line_width=0,
            name="Reverse strand region",
            showlegend=(i == 0),
        )

    fig.update_layout(
        dict(
            title_text=f"SNPs visualization ({snps_path.parts[-2]})",
            legend=dict(
                orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1
            ),
        )
    )
    fig.update_xaxes(
        dict(
            rangeslider=dict(
                visible=True,
                range=total_alignment_range,
            ),
            title="NT position (In reconstructed sequence)",
            showgrid=True,
            gridcolor="Gray",
            gridwidth=1,
            zeroline=False,
            griddash="solid",
        )
    )
    fig.update_yaxes(
        dict(
            range=(MISMATCH_GRAPH_VALUE - 0.5, UNALIGNED_GRAPH_VALUE + 0.5),
            fixedrange=True,
            title="",
            showgrid=True,
            gridcolor="LightGray",
            gridwidth=1,
            zeroline=False,
            griddash="dash",
            tickvals=[
                UNALIGNED_GRAPH_VALUE,
                MATCH_GRAPH_VALUE,
                GAP_GRAPH_VALUE,
                MISMATCH_GRAPH_VALUE,
            ],
            ticktext=["Unaligned", "Match", "Gap", "Mismatch"],
        )
    )

    return fig


def generate_graphs(
    logger: logging.Logger,
    config: PipelineConfig.VisualizationConfig,
    graphs: tuple[str, ...],
    snps_path: Path,
    one_to_one_coords_path: Path,
    export_path: Path,
) -> VisualizationResults:
    results = VisualizationResults()

    export_path.mkdir(parents=True)

    # Plotly config used across all figures
    plotly_config: dict = dict(
        displaylogo=False,
    )

    color_palette: list[str] = []
    match config.color_palette:
        case "custom":
            color_palette = CUSTOM_PALETTE
        case "okabeito":
            color_palette = OKABEITO_PALETTE

    if "SNPs" in graphs:
        snps_graph_path: Path = export_path / "snps_graph.html"
        snps_fig: Figure | None = generate_snps_graph(
            logger, color_palette, snps_path, one_to_one_coords_path
        )
        if snps_fig is not None:
            snps_fig.write_html(
                file=snps_graph_path,
                config=plotly_config,
                include_plotlyjs=config.standalone_graphs,
            )
            results.snps_graph_path = snps_graph_path
        else:
            logger.error("SNPs graph couldn't be generated")

    return results
