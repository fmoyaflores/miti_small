#!/usr/bin/env python
# # Daisy Algorithm for ImmeDB results
#
# 1. For a new row of '20230614_AMR in TE_Abricate-report' googlesheet, open the ImmeDB results for the strain that corresponds to the row.
# 2. In the ImmeDB result file, search hits with "q.start" value smaller than "START" value of Abricate (Column "C" in '20230614_AMR in TE_Abricate-report'; mobile genetic element region should start prior to antibiotic resistance gene).
# 3. Among those, search hits with "q.end" value larger than "END" value of Abricate (Column "D" in '20230614_AMR in TE_Abricate-report'; mobile genetic element region should end after the antibiotic resistance gene).
# 4. For a hit that satisfies both 2&3, get "subject acc.ver" in ImmeDB results file (ex; NZ_NMTU01000006.1:71946-83411).
# 5. Search the "subject acc.ver" in "Data1_MGE_sequences.fasta" file to figure out the name of the hit (ex;  IMEs459).
#     - If there are multiple hits, any ICE elements take priority to document.
#     - If there's no hit, put 'na' in column O (ImmeDB values) of '20230614_AMR in TE_Abricate-report' googlesheet
# 6. Copy the name of the hit to column O (ImmeDB values) of '20230614_AMR in TE_Abricate-report' googlesheet.
# 7. Copy the corresponding information from the row in the ImmeDB result file (from column B-P) and paste it to column P-AD of '20230614_AMR in TE_Abricate-report' googlesheet.

# USAGE: python integrate_Abricate_ImmeDB.py --abricate <abricate result file> --immedb_annotations <immedb annotations file> --output <output file> --project_folder <project folder>
# Example: python integrate_Abricate_ImmeDB.py --abricate abricate_output.csv --output abricate_immedb_output.csv

import argparse
import logging
from collections import Counter

import pandas as pd
from cloudpathlib import AnyPath

logger = logging.getLogger(__name__)


def read_immedb_result(immedb_file, columns):
    return pd.read_table(
        immedb_file,
        comment="#",
        header=None,
        names=columns,
        dtype={"query": str, "subject": str},
    )


def daisy_filter(
    genome_name: str,
    seq_name: str,
    start: int,
    end: int,
    immedb_annotations: dict,
    immedb_columns: list,
    annotation_columns: list,
    project_folder: str,
):
    # s3://genomics-workflow-core/Results/Blast/MITI-MCB/SH0001342-00095/immeDB/SH0001342-00095.blastn.tsv
    immedb_file = AnyPath(
        f"{project_folder}/{genome_name}/immeDB/{genome_name}.blastn.tsv"
    )
    logger.info(
        f"Applying filters to immedb results for '{genome_name}' from '{immedb_file}'"
    )
    df = (
        read_immedb_result(immedb_file, immedb_columns)
        .query("(query == @seq_name) and (q_start <= @start) and (q_end >= @end)")
        .assign(ImmeDB_values=lambda x: x["subject"].map(immedb_annotations))
    )

    if df.empty:
        return [None] * (len(immedb_columns) + len(annotation_columns))

    # add freq of also found categories.
    uniq_annot = Counter(df["ImmeDB_values"].to_list())
    df_top_hit = df.iloc[0].copy()
    df_top_hit["All_immeDB_annotations"] = uniq_annot

    return df_top_hit


def read_abricate_result(abricate_file):
    return (
        pd.read_csv(
            abricate_file,
            sep=",",
            header=0,
            index_col=False,
            dtype=str,
            usecols=set(range(14)),
        )
        .rename(columns={"#FILE": "FILE"})
        .sort_values(by=["FILE"])
    )


def read_immedb_annotations(immedb_annotations_file):
    df = pd.read_csv(
        immedb_annotations_file,
        sep=",",
        header=0,
        index_col=False,
        dtype=str,
        names=["accession", "annotation"],
    )

    return dict(zip(df["accession"], df["annotation"]))


# write an argument parser
def add_arguments():
    parser = argparse.ArgumentParser(
        description="Integrate Abricate and ImmeDB results",
        epilog="""
        Example:
            python integrate_Abricate_ImmeDB.py --abricate abricate_output.csv --output abricate_immedb_output.csv
        """,
    )
    # add usage and example

    parser.add_argument(
        "--abricate",
        type=str,
        help="abricate result file",
    )
    parser.add_argument(
        "--immedb_annotations",
        type=str,
        help="immedb annotations file",
        default="s3://genomics-workflow-core/scratch/sunitj/daisy/immeDB/Data1_MGE_sequences.annotations.csv",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="output file",
    )
    parser.add_argument(
        "--project_folder",
        type=str,
        help="project folder",
        default="s3://genomics-workflow-core/Results/Blast/MITI-MCB",
    )

    return parser.parse_args()


def main():
    args = add_arguments()
    abricate_file = AnyPath(args.abricate)
    immedb_annotations_file = AnyPath(args.immedb_annotations)
    output_file = AnyPath(args.output)
    project_folder = AnyPath(args.project_folder)

    logger.info(f"Reading abricate results file: {abricate_file}")
    abricate_df = read_abricate_result(abricate_file)

    logger.info(f"Reading immedb annotations file: {immedb_annotations_file}")
    immedb_annotations = read_immedb_annotations(immedb_annotations_file)

    immedb_columns = [
        "query",
        "subject",
        "perc_identity",
        "alignment_length",
        "mismatches",
        "gap_opens",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "evalue",
        "bit_score",
        "query_length",
        "subject_length",
        "perc_query_coverage_per_subject",
        "subject sci names",
    ]
    annotation_columns = [
        "ImmeDB_values",
        "All_immeDB_annotations",
    ]

    result_columns = immedb_columns + annotation_columns
    abricate_df[result_columns] = abricate_df.apply(
        lambda x: read_immedb_result(
            AnyPath(x.FILE).stem,
            x.SEQUENCE,
            int(x.START),
            int(x.END),
            immedb_annotations,
            immedb_columns,
            annotation_columns,
            project_folder,
        ),
        axis=1,
        result_type="expand",
    )

    abricate_df.to_csv(
        output_file,
        index=False,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    main()
