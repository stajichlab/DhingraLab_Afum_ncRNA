#!/usr/bin/env python3

import argparse
import sys
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="Process some input.")
    parser.add_argument('-i', '--input', nargs='*', type=str, help='Input files or stdin if not provided')
    parser.add_argument('-o', '--output', type=str, help='Output file or stdout if not provided')
    return parser.parse_args()

def process_blast(input_data,out):
    # Placeholder for processing logic
    incsv = csv.reader(input_data.splitlines(), delimiter='\t')
    outcsv = csv.writer(out, delimiter='\t')
    for row in incsv:
        if row[0].startswith("#"):
            continue
        gene = row[0]
        chrom = row[1]
        gene_start = row[6]
        gene_end = row[7]
        strand = '+'
        if int(gene_start) > int(gene_end):
            strand = '-'
        chrom_start = row[8]
        chrom_end = row[9]
        outcsv.writerow([gene,chrom,chrom_start,chrom_end,strand])

def main():
    args = parse_args()
    
    if args.input:        
        input_data = ""
        for input_file in args.input:
            with open(input_file, 'r') as infile:
                input_data += infile.read()
    else:
        input_data = sys.stdin.read()

    outfh = sys.stdout
    if args.output:
        with open(args.output, 'w') as outfile:
            outfh = outfile

    result = process_blast(input_data, outfh)
    
if __name__ == "__main__":
    main()