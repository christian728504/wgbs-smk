import polars as pl
import os

wildcard_constraints:
    barcode="[A-Za-z0-9]+",
    contig="(chr[0-9XY]+|Pool)"

if not os.path.exists("results/logfiles"):
    os.makedirs("results/logfiles")

tmpdir = config["tmpdir"]
METADATA = pl.read_csv(config["metadata"], separator="\t", comment_prefix="#")

barcodes = METADATA["Barcode"].to_list()
datasets = METADATA["Dataset"].to_list()

experiments = METADATA.select("Dataset").with_columns(
    pl.col("Dataset").str.split("_").list.get(0).alias("Experiment")
)["Experiment"].unique().to_list()

read_1_files = METADATA["File1"].to_list()
read_2_files = METADATA["File2"].to_list()
CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "Pool"]

dataset_lookup = dict(zip(barcodes, datasets))
READ_LOOKUP = dict(zip(barcodes, zip(read_1_files, read_2_files)))

def get_r1(wildcards):
    r1 = READ_LOOKUP.get(wildcards.barcode)[0]
    return r1

def get_r2(wildcards):
    r2 = READ_LOOKUP.get(wildcards.barcode)[1]
    return r2

mapping_patterns = [
    "results/mapping/{barcode}/{barcode}.bam",
    "results/mapping/{barcode}/{barcode}.bam.md5", 
    "results/mapping/{barcode}/{barcode}.json",
    "results/mapping/{barcode}/.continue"
]

chrom_ctg_bed_pattern = "results/calls/{barcode}/{barcode}_{contig}_ctgs.bed"

def get_chrom_bcfs(wildcards):
     return expand("results/calls/{barcode}/{barcode}_{contig}.bcf", barcode=wildcards.barcode, contig=CONTIGS)

calling_patterns = [
    "results/calls/{barcode}/{barcode}_{contig}.bcf",
    "results/calls/{barcode}/{barcode}_{contig}.json"
]

def make_bedmethyl_pearson_correlation(wildcards):
    grouped = (
        METADATA
        .with_columns(
            pl.col("Dataset").str.split("_").list.get(0).alias("Experiment")
        )
        .group_by("Experiment")
        .agg("Barcode")
        .filter(pl.col("Experiment") == wildcards.experiment)
        .select("Barcode")
        .item()
    )
    inputs = [f"results/extract/{barcode}/{barcode}_cpg.bed.gz" for barcode in grouped]
    return inputs
