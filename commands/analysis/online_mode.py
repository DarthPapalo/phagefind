import logging
import re
import time
from typing import TypedDict, cast

import requests

from _commons import PROGRAM_NAME

BLAST_BASE_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
PROKARIOTE_DB = "nt_prok"
# Title:Prokaryota (bacteria and archaea) nt
# Description:The Prokaryote nucleotide collection consists of GenBank+EMBL+DDBJ+PDB+RefSeq sequences,
# but excludes EST, STS, GSS, WGS, TSA, patent sequences as well as phase 0, 1, and 2 HTGS sequences
# and sequences longer than 100Mb. The database is non-redundant.
# Identical sequences have been merged into one entry,while preserving the accession, GI, title and
# taxonomy information for each entry.
# Molecule Type:mixed DNA

ESUMMARY_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

POLLING_TIME = 66


class BlastSearchResult(TypedDict):
    identified_phage_title: str
    new_identified_targets: dict[
        str, tuple[str, str, float, int]
    ]  # ACC : (Organism name, strain, lowest evalue, highest score)


def make_blast_request(
    logger: logging.Logger, phage_id: str, email: str
) -> BlastSearchResult | None:
    """
    Returns `None` if any error occurred during the request.
    """
    put_request = f"?CMD=Put&EMAIL={email}&TOOL={PROGRAM_NAME}&PROGRAM=blastn&MEGABLAST=on&DATABASE={PROKARIOTE_DB}&QUERY={phage_id}"
    put_result: requests.Response = requests.post(BLAST_BASE_URL + put_request)

    if put_result.status_code != 200:
        logger.error(
            f"Error making POST HTTP request for online-mode analysis. Status code [red]{put_result.status_code}[/red]"
        )
        return None

    # Parse Job RID from response
    rid_match: re.Match[str] | None = re.search(
        r"^    RID = (.*$)", put_result.text, re.MULTILINE
    )
    if rid_match is None:
        logger.error("Error parsing BLAST Job RID from request response")
        return None
    rid: str = rid_match.group(1)

    # Parse estimated time for completion
    rtoe_match: re.Match[str] | None = re.search(
        r"^    RTOE = (.*$)", put_result.text, re.MULTILINE
    )
    if rtoe_match is None:
        logger.error("Error parsing BLAST Job RTOE from request response")
        return None
    rtoe: str = rtoe_match.group(1)

    logger.info(f"Successful BLAST request. Estimated completion time: {rtoe}")

    # Wait for request completion
    time.sleep(int(rtoe) + 10)

    # Poll for completion
    total_poll_time: int = int(rtoe) + 10
    while True:
        get_request = f"?CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular&RID={rid}"
        get_result: requests.Response
        for _attempt in range(3):
            try:
                get_result: requests.Response = requests.get(
                    BLAST_BASE_URL + get_request
                )
                break
            except requests.exceptions.ChunkedEncodingError:
                time.sleep(1)
        else:
            logger.error(f"Failed to retrieve BLAST results after {_attempt} attempts")
            return None

        if get_result.status_code != 200:
            logger.error(
                f"Error making GET HTTP request for online-mode analysis. Status code {get_result.status_code}"
            )
            return None

        status_match: re.Match[str] | None = re.search(
            r"\s+Status=(\w+)", get_result.text, re.MULTILINE
        )

        if status_match is None:
            logger.error(
                f"Error parsing BLAST Job status from request response (RID: {rid})"
            )
            return None
        status: str = status_match.group(1)
        match status:
            case "WAITING":
                logger.debug(f"polling... (RID: {rid}, total time: {total_poll_time})")
                time.sleep(POLLING_TIME)
                total_poll_time += POLLING_TIME
                continue
            case "FAILED":
                logger.error(f"BLAST search for phage {phage_id} failed (RID: {rid})")
                return None
            case "UNKNOWN":
                logger.error(f"BLAST search for phage {phage_id} expired (RID: {rid})")
                return None
            case "READY":
                query_name: str = ""
                best_organism_evalue_and_bitscore: dict[str, tuple[float, int]] = {}
                acc_to_organisms_strains: dict[str, tuple[str, str]] = {}
                number_hits: int = -1

                # Parse blast result
                for line in get_result.text.splitlines():
                    # 1. Parse query name
                    if len(query_name) < 1:
                        query_name_match: re.Match[str] | None = re.search(
                            r"^# Query=\s*(.+)$", line
                        )
                        if query_name_match is None:
                            continue
                        else:
                            query_name = query_name_match.group(1)

                    # 2. Parse number of hits
                    if number_hits < 0:
                        # "No significant similarity found.""
                        number_hits_match: re.Match[str] | None = re.search(
                            r"^# (\d+) hits found", line
                        )
                        if number_hits_match is None:
                            continue
                        else:
                            logger.info(
                                f"BLAST search for phage {phage_id} completed, {number_hits} hits, retrieving results..."
                            )
                            number_hits = int(number_hits_match.group(1))

                    # 3. Parse hits entries
                    if number_hits > 0:
                        (
                            query_acc,
                            subject_acc,
                            identity,
                            alignment_length,
                            mismatches,
                            gaps,
                            query_start,
                            query_end,
                            subject_start,
                            subject_end,
                            evalue,
                            bit_score,
                        ) = line.split("\t")

                        converted_evalue = float(evalue)
                        converted_bitscore = int(bit_score)
                        current_values: tuple[float, int] | None = (
                            best_organism_evalue_and_bitscore.get(subject_acc)
                        )
                        if current_values is not None:
                            # Prefer lower evalues even if it means lower bit scores
                            if converted_evalue > current_values[0]:
                                continue
                            elif (
                                converted_evalue == current_values[0]
                                and converted_bitscore < current_values[1]
                            ):
                                continue
                        best_organism_evalue_and_bitscore[subject_acc] = (
                            converted_evalue,
                            converted_bitscore,
                        )

                    # 4. Get organisms names and strains for accessions
                    summary_request = f"?db=nuccore&id={','.join(best_organism_evalue_and_bitscore.keys())}&retmode=json"
                    summary_result: requests.Response = requests.get(
                        ESUMMARY_BASE_URL + summary_request
                    )

                    accesions_summaries = cast(
                        dict[str, dict], summary_result.json()["result"]
                    )
                    accesions_summaries.pop("uids")

                    for uid in accesions_summaries:
                        acc: str = accesions_summaries[uid]["accessionversion"]
                        acc_to_organisms_strains[acc] = (
                            accesions_summaries[uid]["organism"],
                            accesions_summaries[uid]["strain"],
                        )

                    return BlastSearchResult(
                        identified_phage_title=query_name,
                        new_identified_targets={
                            acc: (
                                *acc_to_organisms_strains[acc],
                                *best_organism_evalue_and_bitscore[acc],
                            )
                            for acc in best_organism_evalue_and_bitscore
                        },
                    )

            case _:
                logger.error(
                    f"Unknown BLAST request status code received: {status}. Aborting online mode analysis"
                )
                return None
