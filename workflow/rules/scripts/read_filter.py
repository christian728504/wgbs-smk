import sys
import pysam
import argparse
from typing import Dict, Optional

# @SQ     SN:chr2 LN:242193529    M5:f98db672eb0993dcfdabafe2a882905c     UR:file:///zata/zippy/ramirezc/gembs-smk/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz        AS:GRCh38       SP:Homo_sapiens
# @SQ     SN:chr3 LN:198295559    M5:76635a41ea913a405ded820447d067b0     UR:file:///zata/zippy/ramirezc/gembs-smk/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz        AS:GRCh38       SP:Homo_sapiens

def read_contig_file(filename: str) -> Dict[str, Dict[str, str]]:
    """Read contig metadata file and return as nested dict."""
    contigs = {}
    with open(filename, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if not fields:
                continue
            
            name = fields[0]
            tags = {}
            
            for field in fields[1:]:
                if field.startswith(('LN:', 'M5:', 'AS:', 'SP:')):
                    tag_type = field[:2]
                    tag_value = field[3:]
                    tags[tag_type] = tag_value
            
            contigs[name] = tags
    
    return contigs

def clean_read_name(read_name: str) -> str:
    """Clean read name by replacing invalid characters."""
    cleaned = []
    for char in read_name:
        if char == '@':
            cleaned.append('_')
        elif 33 <= ord(char) <= 126:
            cleaned.append(char)
        else:
            cleaned.append('_')
    
    return ''.join(cleaned)

def update_header(header: pysam.AlignmentHeader, 
                  contig_metadata: Optional[Dict[str, Dict[str, str]]]) -> pysam.AlignmentHeader:
    """Update @SQ header lines with contig metadata."""
    if not contig_metadata:
        return header
    
    new_header = header.to_dict()
    
    if 'SQ' in new_header:
        for sq_record in new_header['SQ']:
            seq_name = sq_record.get('SN')
            if seq_name and seq_name in contig_metadata:
                for tag, value in contig_metadata[seq_name].items():
                    sq_record[tag] = value
    
    return pysam.AlignmentHeader.from_dict(new_header)

def main():
    parser = argparse.ArgumentParser(description='Filter SAM/BAM read names and update headers')
    parser.add_argument('contig_file', nargs='?', help='Optional contig metadata file')
    parser.add_argument('-i', '--input', default='-', help='Input SAM/BAM file (default: stdin)')
    parser.add_argument('-o', '--output', default='-', help='Output SAM/BAM file (default: stdout)')
    
    args = parser.parse_args()
    
    contig_metadata = None
    if args.contig_file:
        contig_metadata = read_contig_file(args.contig_file)
    
    input_file = sys.stdin if args.input == '-' else args.input
    with pysam.AlignmentFile(input_file, 'r') as infile:
        
        updated_header = update_header(infile.header, contig_metadata)
        
        output_file = sys.stdout if args.output == '-' else args.output
        output_mode = 'w' if args.output == '-' else 'wb'
        
        with pysam.AlignmentFile(output_file, output_mode, header=updated_header) as outfile:
            
            for read in infile:
                read.query_name = clean_read_name(read.query_name)
                outfile.write(read)

if __name__ == '__main__':
    main()
