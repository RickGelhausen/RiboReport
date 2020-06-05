#!/usr/bin/env python
#./filter.py --in_gff_filepath escherichia_coli/annotation.gtf --in_locus_tag_filepath escherichia_coli/ecoli_massspec.tsv > ecoli_masspec.gff

import argparse
import os
import csv
import re

def get_tag_set(input_csv_filepath):
    #print("Getting tags:")
    tag_set = set()
    with open(input_csv_filepath, newline='\n') as csvfile:
        for tag in csv.reader(csvfile):
            #print(tag)
            if not len(tag) == 0:
                tag_set.add(tag[0])
        #print("Computed tag set")
    return tag_set

def filter_gff(input_gff_filepath, tag_set, filter_tag, feature_type):
    #entries = []
    with open(input_gff_filepath, newline='\n') as csvfile:
        gffreader = csv.reader(csvfile, delimiter='\t')
        for entry in gffreader:
            if (not entry[0].startswith('#')):
                gfftype = entry[2]
                if gfftype == feature_type:
                    seqid = entry[0]
                    source = entry[1]
                    gfftype = entry[2]
                    start = entry[3]
                    end = entry[4]
                    score = entry[5]
                    strand = entry[6]
                    phase = entry[7]
                    attributes = entry[8]
                    if filter_tag == "gene":
                        gene = re.search(r"gene=\w+",attributes)
                        gene_id =  re.sub('^gene=', '', gene.group(0))
                        #print(gene_id)
                        if gene_id in tag_set:
                            print("\t".join(entry))
                    if filter_tag == "locus_tag":
                        locus = re.search(r"locus_tag=\w+",attributes)
                        locus_tag =  re.sub('^locus_tag=', '', locus.group(0))
                        #print(gene_id)
                        if locus_tag in tag_set:
                            print("\t".join(entry))


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Filter GFF with list of locus tags')
    parser.add_argument("--in_gff_filepath", help='Input wig filepath', required=True)
    parser.add_argument("--in_locus_tag_filepath", help='Input wig filepath', required=True)
    parser.add_argument("--feature_type", help="Feature type to extract, e.g. gene or CDS", required=True)
    parser.add_argument("--filter_tag", help='Tag to filter for, either locus_tag or gene', required=True)
    args = parser.parse_args()
    locus_set = get_tag_set(args.in_locus_tag_filepath)
    filter_gff(args.in_gff_filepath, locus_set, args.filter_tag, args.feature_type)

if __name__ == '__main__':
    main()
