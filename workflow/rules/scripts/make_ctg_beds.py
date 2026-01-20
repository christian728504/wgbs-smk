import polars as pl
import os
import argparse

def parse_contig_sizes(contig_file: str):
    """Parse the contig sizes file into a pandas DataFrame"""
    df = pl.read_csv(contig_file, separator='\t', has_header=False, new_columns=['chrom', 'end'])
    df = df.with_columns(
        pl.lit(0).alias('start')
    ).select("chrom", "start", "end")
    return df

def create_bed_files(df: pl.DataFrame, output_dir: str, barcode: str):
    """Create BED files for each chromosome and Pool@1"""
    os.makedirs(output_dir, exist_ok=True)
    main_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    
    for chrom in main_chroms:
        chrom_data = df.filter(pl.col('chrom') == chrom)
        if not chrom_data.is_empty():
            bed_file = os.path.join(output_dir, f'{barcode}_{chrom}_ctgs.bed')
            chrom_data.write_csv(bed_file, include_header=False, separator="\t")
            print(f"Created: {bed_file}")
    
    alt_contigs = df.filter(~(pl.col('chrom').is_in(main_chroms)))
    if not alt_contigs.is_empty():
        pool_file = os.path.join(output_dir, f'{barcode}_Pool_ctgs.bed')
        alt_contigs.write_csv(pool_file, include_header=False, separator="\t")
        print(f"Created: {pool_file}")
        print(f"Pool contains {len(alt_contigs)} alternative contigs")

def main():
    parser = argparse.ArgumentParser(description='Generate contig BED files from contig sizes')
    parser.add_argument('--chromsizes', help='Path to contig sizes file')
    parser.add_argument('--output_dir', help='Output directory for BED files')
    parser.add_argument('--barcode', help='Sample barcode/name')
    
    args = parser.parse_args()
    
    print(f"Reading chromsizes from: {args.chromsizes}")
    df = parse_contig_sizes(args.chromsizes)
    print(f"Found {df.shape[0]} contigs")
    
    create_bed_files(df, args.output_dir, args.barcode)
    print("Done!")

if __name__ == '__main__':
    main()
