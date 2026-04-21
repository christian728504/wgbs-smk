# WGBS Snakemake Pipeline

## Overview

Here we present a [snakemake](https://github.com/snakemake/snakemake.git) pipeline for whole genome bisulfite sequencing (WGBS) data processing, including read mapping, methylation calling, and signal extraction across multiple samples.

> [!NOTE]
> This pipeline uses [gemBS-rs](https://github.com/heathsc/gemBS-rs) and [GEM3](https://github.com/smarco/gem3-mapper) for bisulfite-aware alignment and methylation calling.

## Example Usage

### Prerequisites

To run this pipeline, you will need paired-end WGBS FASTQ files and a reference genome (e.g., GRCh38). You will also need a tab-delimited metadata file describing your samples. See [Metadata File](#metadata-file) for details.

### Dependencies

This pipeline requires the following dependencies:

- [`conda` (preferably `mamba`)](https://github.com/conda-forge/miniforge)
- [`singularity`](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
- (Optional) [`slurm`](https://slurm.schedmd.com/quickstart.html)

To start, first setup the runtime environment:

```bash
mamba env create --name=wgbs-smk --file=environments/wgbs-smk.yaml
mamba activate wgbs-smk
```

### Metadata File

The pipeline requires a tab-delimited metadata file with a header row. At minimum, the file must contain the following columns:

- `Barcode`: A unique sample identifier.
- `Dataset`: The dataset or sequencing run name.
- `File1`: The path to the forward (R1) FASTQ file.
- `File2`: The path to the reverse (R2) FASTQ file.

Lines beginning with `#` are treated as comments and ignored.

### Example Configfile

For a complete example, see [config/config.yml](config/config.yml). From this file you'll want to pay attention to the following sections:

```yaml
workdir: "/data/zusers/ramirezc/wgbs-smk"
metadata: "jobs/2026-02-14/metadata.tsv"
reference: "jobs/2026-02-14/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
assembly: "GRCh38"
species: "Homo_sapiens"
project: "MOHD"
tmpdir: "/tmp"
make_index: True
```

- `workdir`: The working directory for the pipeline. All relative paths in the config are resolved relative to this directory.
- `metadata`: Path to the tab-delimited metadata file describing your samples.
- `reference`: Path to the reference genome FASTA file. Gzipped FASTAs are supported and will be decompressed automatically.
- `assembly`: (Optional) Assembly name, used in the SAM sequence dictionary header.
- `species`: (Optional) Species name (no whitespace), used in the SAM sequence dictionary header.
- `project`: Project name used for output file naming.
- `tmpdir`: Temporary directory for intermediate files. In some compute environments, `/tmp` is mounted with `noexec`. In this case, you may need to use a different location. Note that the SLURM profile uses `shadow-prefix: /tmp` for shadow rules.
- `make_index`: Whether to build the reference index. If you have pre-built indexes, set this to `False`.

```yaml
index:
  sampling_rate: 4
fastq_split: 2
```

- `index.sampling_rate`: The text sampling rate for GEM indexing. Lower values produce larger but faster indexes.
- `fastq_split`: The number of parts to split each FASTQ pair into for parallel mapping. Each part is mapped independently and merged afterwards.

```yaml
mapping:
  non_stranded: false
  remove_individual_bams: true
```

- `mapping.non_stranded`: Whether to use non-stranded mapping mode.
- `mapping.remove_individual_bams`: Whether to remove individual per-split BAM files after merging.

```yaml
calling:
  output_type: "b"
  mapq_threshold: 10
  qual_threshold: 13
  reference_bias: 2
  left_trim: 5
  right_trim: 0
  keep_improper_pairs: false
  keep_duplicates: false
  haploid: false
  conversion: 0.01,0.05
  remove_individual_bcfs: true
  contig_pool_limit: 25000000
```

- `calling.output_type`: Output format for `bs_call` (e.g., `b` for BCF, `v` for VCF).
- `calling.mapq_threshold`: Minimum mapping quality for a read to be used in calling.
- `calling.qual_threshold`: Minimum base quality threshold.
- `calling.reference_bias`: Reference bias parameter for `bs_call`.
- `calling.left_trim` / `calling.right_trim`: Number of bases to trim from the left and right ends of reads during calling.
- `calling.conversion`: Bisulfite conversion rate parameters (under-conversion, over-conversion), as a comma-delimited string.
- `calling.contig_pool_limit`: Maximum size (in bp) for contigs to be pooled together in the `Pool` contig group.

```yaml
extract:
  strand_specific: true
  bw_strand_specific: true
  phred_threshold: 10
```

- `extract.strand_specific`: Whether to extract methylation in strand-specific mode.
- `extract.bw_strand_specific`: Whether to produce strand-specific BigWig files.
- `extract.phred_threshold`: Minimum Phred quality threshold for extraction.

### Defining Resources

All resource definitions for rules are located in [profile/slurm/config.yaml](profile/slurm/config.yaml).

```yaml
  mapping:
    runtime: 7200m
    threads: 26
    cpus_per_task: 26
    slurm_partition: 5days
    constraint: cascadelake
    slurm_extra: "--exclude=z[1024,1062] --cpu-freq=High-High:Performance"
    mem: 250000MB
```

In this instance, we are requesting a job on a SLURM cluster with a runtime of 7200 minutes, a specific constraint (i.e. cascadelake), a specific partition (i.e. 5days), 26 threads, 26 CPUs per task, and 250000MB of memory.

These resource definitions are specific to the Weng Lab's SLURM cluster. If you are running the pipeline on your own cluster, you'll need to adjust these values accordingly. We recommend consulting the [snakemake-executor-plugin-slurm documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) for further information.

### Running the Pipeline

With the hard part out of the way, you can now run the pipeline:

```bash
snakemake --workflow-profile profile/slurm
```

> [!NOTE]
> It is recommended that you run the pipeline in dry run mode first (add the `-n` flag). Also note that `snakemake` must be ran in the root of this repository.

This will run the pipeline and produce all outputs to a directory `results`.

### What if I cannot run on a SLURM cluster?

If you cannot run the pipeline on a SLURM cluster, you can run it locally using the default profile:

```bash
snakemake --workflow-profile profile/default
```

The default profile is configured for local execution. Resource definitions for the local profile are located in [profile/default/config.yaml](profile/default/config.yaml). Note that this pipeline is resource-intensive and local execution may not be feasible for large datasets.

## Outputs

The pipeline produces results organized under the `results/` directory:

| Directory | Key Outputs |
|---|---|
| `results/indexes/` | Reference FASTA, FASTA index, chromosome sizes TSV, SAM sequence dictionary, GEM bisulfite index (`.gem`) |
| `results/split/` | Per-barcode split paired-end FASTQs for parallel mapping |
| `results/mapping/` | Per-split BAMs, merged sorted BAM, mapping statistics JSON |
| `results/postprocess/` | Filtered and sorted BAM, CSI index, MD5 checksum |
| `results/get_coverage/` | Strand-specific coverage bedGraphs (`.bg`) and BigWig files (`.bw`) |
| `results/calls/` | Per-contig BCFs and JSON reports, merged whole-genome BCF with CSI index and MD5 checksum |
| `results/extract/` | CpG and non-CpG methylation text files (`.txt.gz`), bedMethyl files for CpG/CHG/CHH contexts (`.bed.gz`), strand-specific BigWig files |
| `results/extract_signal/` | Strand-specific BigWig files for CpG, CHG, and CHH contexts (e.g., `{barcode}_cpg_pos.bw`, `{barcode}_chh_neg.bw`) |
| `results/get_methylc/` | Combined methylC BED file (`.bed.gz`) with tabix index, containing CpG, CHG, and CHH sites |

## Pipeline Steps

1. **Indexing** (`indexing.smk`): Prepares the reference genome by creating a FASTA index, chromosome sizes file, SAM sequence dictionary, and GEM bisulfite index.
2. **Split** (`mapping.smk`): Splits paired-end FASTQ files into parts for parallel mapping using `seqkit split2`.
3. **Mapping** (`mapping.smk`): Aligns bisulfite-converted reads to the reference genome using `gem-mapper` with bisulfite mode.
4. **Merge** (`mapping.smk`): Merges per-split BAM files and sorts the combined alignment.
5. **Postprocess** (`mapping.smk`): Filters reads, recalculates coverage statistics, and produces final BAM files with checksums.
6. **Coverage** (`mapping.smk`): Generates strand-specific coverage bedGraphs and BigWig files using `bedtools genomecov` and `bedGraphToBigWig`.
7. **Generate Contig BEDs** (`call.smk`): Creates per-chromosome and pooled contig BED files for parallelized methylation calling.
8. **Calling** (`call.smk`): Calls methylation at each CpG/CHG/CHH site per contig using `bs_call`.
9. **Merge Chroms** (`call.smk`): Concatenates per-contig BCF files into a single whole-genome BCF per sample.
10. **Extract** (`extract.smk`): Extracts methylation information from BCF files using `mextr`, producing bedMethyl files, compressed CpG/non-CpG text files, and strand-specific BigWig files.
11. **Extract Signal** (`extract.smk`): Converts bedMethyl files to strand-specific BigWig files for CpG, CHG, and CHH contexts.
12. **Get methylC** (`extract.smk`): Combines CpG, CHG, and CHH bedMethyl data into a single sorted, bgzipped, and tabix-indexed methylC BED file per sample.

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

> Marco-Sola, S., Sammeth, M., Guigó, R., & Ribeca, P. _The GEM mapper: fast, accurate and versatile alignment by filtration_. Nature Methods, 9(12), 1185–1188, **2012**. https://doi.org/10.1038/nmeth.2221.

> https://github.com/heathsc/gemBS-rs

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
