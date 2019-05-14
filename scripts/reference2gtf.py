#!/usr/bin/env python
import csv
#import pandas as pd
#import re
import argparse
#import numpy as np
#import os

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Converts reference to gtf.')
    parser.add_argument("--input_tsv_filepath", help='Path to write tsv output')
    parser.add_argument("--output_gtf_filepath", help='Path to write gtf output')
    args = parser.parse_args()
    tsvin =  open(args.input_tsv_filepath,'r', encoding='utf-8')
    gtfout = open(args.output_gtf_filepath, 'w', encoding='utf-8')
    tsvlines = csv.reader(tsvin, delimiter='\t')
    output=""
    for row in tsvlines:
         seqid = "1"
         start_coordinate = row[0]
         end_coordinate = row[1]
         strand = row[2]
         orf_id = row[3]
         locus_tag= row[5]
         gene_name= row[4]
         transcribed= row[6]
         translated= row[7]
         attributes = "gene_id \"" + gene_name + "\";"
         entry = seqid + "\t" + "annotation" + "\t" + "CDS" + "\t" + str(start_coordinate)  + "\t" + str(end_coordinate) + "\t" + "." + "\t" + strand + "\t" + "." + "\t" + attributes + "\n"
         output+=entry
    gtfout.write(output)

if __name__ == '__main__':
    main()
