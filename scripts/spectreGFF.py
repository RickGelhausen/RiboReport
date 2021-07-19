#!/usr/bin/env python
"""
Script to transform SPECtre output to prediction gtf format in order to run the benchmark evaluation scripts on it.
"""

import pandas as pd
import argparse
import csv
import collections

def spectre_to_gtf(spectre_input, condition):
    """
    Transform spectre output into gtf
    """

    input_df = pd.read_csv(spectre_input, sep="\t", comment="#")

    n_tuple = collections.namedtuple("Pandas", ["chrom", "source", "type", "start", "stop", "score", "strand", "phase", "attribute"])

    rows = []
    for row in input_df.itertuples(index=True, name="Pandas"):
        chromosome = getattr(row, "chr")
        strand = getattr(row, "strand")
        gene_type = getattr(row, "gene_type")
        ribo_fpkm = getattr(row, "ribo_fpkm")
        coordinates_cds = getattr(row, "coordinates_CDS")

        if gene_type != "protein_coding":
            continue

        start, stop = int(coordinates_cds.split("-")[0]), int(coordinates_cds.split("-")[-1])

        identifier = "%s:%s-%s:%s" % (chromosome, start, stop, strand)

        #result = [chromosome, "spectre", "CDS", start+1, stop, ribo_fpkm, strand, ".", "gene_id \"%s\"; method \"%s\"; condition \"%s\"" % (identifier, "spectre", condition)]
        result = [chromosome, "spectre", "CDS", start+1, stop, ribo_fpkm, strand, ".", "ID=%s;Method=%s;Condition=%s" % (identifier, "spectre", condition)]    
        rows.append(n_tuple(*result))

    return pd.DataFrame.from_records(rows, columns=["chrom", "source", "type", "start", "stop", "score", "strand", "phase", "attribute"])

def main():
    parser = argparse.ArgumentParser(description="Convert SPECtre output to gtf format")
    parser.add_argument("-i", "--input_txt", action="store", dest="spectre_input", required=True, help= "input file ( output result file created by spectre)")
    parser.add_argument("-c", "--condition", action="store", dest="condition", required=True, help= "condition of current sample")
    parser.add_argument("-o", "--out_gtf", action="store", dest="out_gtf", required=True, help= "output gtf")

    args = parser.parse_args()

    gtf_df = spectre_to_gtf(args.spectre_input, args.condition)

    gtf_df.to_csv(args.out_gtf, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()

