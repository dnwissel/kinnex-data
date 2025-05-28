import argparse
from typing import List

import pandas as pd
from typeguard import typechecked


@typechecked
def aggregate_counts(
    sample_paths: List[str],
    sample_names: List[str],
    gtf_annotation_path: str,
    output_path_transcript: str,
    output_path_gene: str,
    transcript_id_col_ix: int,
    count_id_col_ix: int,
) -> int:
    df_dict = {}
    for sample_path, sample_name in zip(sample_paths, sample_names):
        df = pd.read_csv(sample_path, sep="\t")
        if sample_name == sample_names[0]:
            df_dict["transcript_id"] = df.iloc[
                :, (transcript_id_col_ix - 1)
            ].values.tolist()
        df_dict[sample_name] = df.iloc[:, (count_id_col_ix - 1)].values.tolist()
    df = pd.DataFrame(df_dict)
    gtf_df = pd.read_csv(gtf_annotation_path, sep="\t", header=None)
    gtf_df.columns = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    gtf_df["transcript_id"] = gtf_df["attribute"].str.extract(
        r'transcript_id "([^"]+)"'
    )
    gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'gene_id "([^"]+)"')
    gene_transcript_df = gtf_df[["gene_id", "transcript_id"]].copy(deep=True)
    del gtf_df
    df = df.merge(
        gene_transcript_df,
        how="inner",
        left_on="transcript_id",
        right_on="transcript_id",
    ).drop_duplicates()

    df = df[["gene_id", "transcript_id"] + sample_names]

    df = df.sort_values(["gene_id", "transcript_id"])
    df.to_csv(path_or_buf=output_path_transcript, sep="\t", header=True, index=False)
    df.groupby("gene_id").aggregate("sum").reset_index().drop(
        columns=["transcript_id"]
    ).to_csv(path_or_buf=output_path_gene, sep="\t", header=True, index=False)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sample_paths",
        nargs="*",
        type=str,
    )
    parser.add_argument(
        "--sample_names",
        nargs="*",
        type=str,
    )
    parser.add_argument("--gtf_annotation_path", type=str)
    parser.add_argument("--output_path_transcript", type=str)
    parser.add_argument("--output_path_gene", type=str)
    parser.add_argument("--transcript_id_col_ix", type=int)
    parser.add_argument("--count_id_col_ix", type=int)

    args = parser.parse_args()
    aggregate_counts(
        sample_paths=args.sample_paths,
        sample_names=args.sample_names,
        gtf_annotation_path=args.gtf_annotation_path,
        output_path_transcript=args.output_path_transcript,
        output_path_gene=args.output_path_gene,
        transcript_id_col_ix=args.transcript_id_col_ix,
        count_id_col_ix=args.count_id_col_ix,
    )
