"""
Calculates Pearson correlation of sites with ≥10X coverage between CpG bedmethyl files.
Inspired by https://github.com/ENCODE-DCC/mirna-seq-pipeline/blob/dev/src/calculate_correlation.py
"""

import argparse
import itertools
import polars as pl
import pandas as pd

DF_COLUMN_NAMES = ("chrom", "start", "end", "coverage", "methylation_percentage")

def load_bedmethyl(path: str) -> tuple[str, pl.DataFrame]:
    basename = path.split("_")[0]
    df = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        skip_rows=1,
        columns=[0, 1, 2, 9, 10],
        new_columns=DF_COLUMN_NAMES
    )
    return (basename, df)

def calculate_pearson(bedmethyl1: pl.DataFrame, bedmethyl2: pl.DataFrame) -> float:
    bedmethyl1 = bedmethyl1.filter(pl.col("coverage") >= 10)
    bedmethyl2 = bedmethyl2.filter(pl.col("coverage") >= 10)
    merged = bedmethyl1.join(bedmethyl2, on=["chrom", "start", "end"])
    pearson_correlation: float = merged.select(pl.corr("methylation_percentage", "methylation_percentage_right", method="pearson")).item()
    return pearson_correlation

def make_pairwise_correlations(bedmethyls: list[tuple[str, pl.DataFrame]]) -> pd.DataFrame:
    product = list(itertools.product(bedmethyls, repeat=2))
    file_names = [basename for (basename, _) in bedmethyls]
    corr_matrix = pd.DataFrame(index=file_names, columns=file_names, dtype=float)
    for (basename1, bedmethyl1), (basename2, bedmethyl2) in product:
        if basename1 == basename2:
            corr_matrix.loc[basename1, basename2] = 1.0
        else:
            corr_matrix.loc[basename1, basename2] = calculate_pearson(bedmethyl1, bedmethyl2)
    return corr_matrix
    
def main(args: argparse.Namespace):
    bedmethyls = [load_bedmethyl(bedmethyl) for bedmethyl in args.bedmethyls]
    corr_matrix = make_pairwise_correlations(bedmethyls)
    corr_matrix.to_csv(args.output_csv, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser() 
    parser.add_argument("--bedmethyls", nargs="+", required=True)
    parser.add_argument("--output-csv", type=str, required=True)
    args = parser.parse_args()
    
    main(args)
