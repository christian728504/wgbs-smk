# Adapted from  https://github.com/ENCODE-DCC/wgbs-pipeline/blob/dev/wgbs_pipeline/calculate_average_coverage.py
"""
Calculates average coverage from a bam file. The formula for this is given by this:
Coverage = (Aligned read counts * read length ) / total genome size

The total genome size is obtained by summing the values in the chromosome sizes. The
aligned read counts is the number of reads in the bam. The read length is obtained from
samtools-stats `RL`, see http://www.htslib.org/doc/samtools-stats.html for details.
"""

import argparse
import tempfile
from typing import Any, Dict, List, Tuple
import sys
import json
import pysam
from pathlib import Path
from collections import OrderedDict
from qc_utils import QCMetric, QCMetricRecord
from qc_utils.parsers import parse_samtools_stats


def get_samtools_stats(bam_path: str, threads: int) -> Dict[str, Any]:
    """
    The buffer must be flushed first before it is read again by the parser, see
    https://stackoverflow.com/questions/46004774/python-namedtemporaryfile-appears-empty-even-after-data-is-written
    """
    samtools_stats_stdout = pysam.stats("--threads", str(threads), bam_path)
    with tempfile.NamedTemporaryFile(mode="w") as samtools_stats_file:
        samtools_stats_file.write(samtools_stats_stdout)
        samtools_stats_file.flush()
        samtools_stats = parse_samtools_stats(samtools_stats_file.name)
    return samtools_stats


def calculate_genome_size(chrom_sizes_path: str) -> int:
    size = 0
    with open(chrom_sizes_path) as f:
        for line in f:
            size += int(line.split()[1])
    return size


def calculate_average_coverage(
    genome_size: int, aligned_read_count: int, read_length: int
) -> Dict[str, float]:
    average_coverage = (aligned_read_count * read_length) / genome_size
    return {"average_coverage": average_coverage}


def make_qc_record(qcs: List[Tuple[str, Dict[str, Any]]]) -> QCMetricRecord:
    qc_record = QCMetricRecord()
    for name, qc in qcs:
        metric = QCMetric(name, qc)
        qc_record.add(metric)
    return qc_record


def main():
    parser = argparse.ArgumentParser(description='Gather additional metrics from the alignment file')
    parser.add_argument('--bamfile', help='Path to bam file')
    parser.add_argument('--chromsizes', help="Path to chromsizes file")
    parser.add_argument('--threads', help='Number of threads (samtools flagstat)')
    parser.add_argument('--gem_mapper_json', help='Where to output metrics')
    args = parser.parse_args()
    
    samtools_stats = get_samtools_stats(args.bamfile, args.threads)
    genome_size = calculate_genome_size(args.chromsizes)
    average_coverage = calculate_average_coverage(
        genome_size=genome_size,
        aligned_read_count=samtools_stats["reads mapped"],
        read_length=samtools_stats["average length"],
    )
    qc_record = make_qc_record(
        [("samtools_stats", samtools_stats), ("average_coverage", average_coverage)]
    )
    
    with open(args.gem_mapper_json, "r") as f:
        existing_data = json.load(f, object_pairs_hook=OrderedDict)
    
    existing_data.update(qc_record.to_ordered_dict())
    
    with open(args.gem_mapper_json, "w") as f:
        json.dump(existing_data, f, indent=2)

if __name__ == "__main__":
    main()
