#!/usr/bin/env python
'''This script takes input files generated by
reparation and creates a new data frame containing specified
information and writes it as gff3 format files.
'''

import pandas as pd
import re
import argparse
import numpy as np
import os
import csv
import collections


# little helper function to create named tuple without having to always state every argument
def createNTuple(args, row):
    nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])
    # txt file content
    ORF_locus = getattr(row, "ORF_locus")
    strand = getattr(row, "strand")
    length = str(getattr(row, "length"))
    start_codon = getattr(row, "start_codon")
    ribo_count = str(getattr(row, "ribo_count"))
    ribo_rpkm = str(getattr(row, "ribo_rpkm"))
    ribo_coverage = str(getattr(row, "ribo_coverage"))
    SD_score = str(getattr(row, "SD_score"))
    SD_pos = str(getattr(row, "SD_pos"))
    prob = str(getattr(row, "prob"))
    ORF_type = getattr(row, "ORF_type")
    Reference = str(getattr(row, "Reference"))
    Distance_from_aTIS = str(getattr(row, "Distance_from_aTIS"))

    # new content
    chromosome, rest = ORF_locus.split(":")
    start, stop = rest.split("-")
    # modify coordinates to include stop codon
    if strand is '+':
       start = str(int(start)) # due to bug in reparation
       stop = str(int(stop) + 3)
    if strand is '-':
       start = str(int(start) - 3)
       stop = str(int(stop)) # due to bug in reparation

    seqName = chromosome
    source = "reparation"
    type = "CDS"
    score = "."
    phase = "."
    attribute = "ID=" + chromosome + ":" + start + "-" + stop + ":" + strand \
              + ";Name=" + chromosome + ":" + start + "-" + stop + ":" + strand \
              + ";ORF_type=" + ORF_type + ";Length=" + length + ";Ribo_count=" + ribo_count \
              + ";Ribo_rpkm=" + ribo_rpkm + ";Ribo_coverage=" + ribo_coverage + ";SD_score=" + SD_score \
              + ";SD_pos=" + SD_pos + ";Prob=" + prob + ";Reference=" + Reference \
              + ";Distance_from_aTIS=" + Distance_from_aTIS + ";Condition=" + args.condition \
              + ";Replicate=" + args.replicate + ";Method=reparation"

    return nTuple(seqName, source, type, start, stop, score, strand, phase, attribute)


def to_gff3(args):
    inputDF = pd.read_csv(args.predictedORFs, sep='\t')

    # extract information from each row and build new dataframe in gff format
    rows = []
    for row in inputDF.itertuples(index=True, name='Pandas'):
        rows.append(createNTuple(args, row))

    return pd.DataFrame.from_records(rows, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"])


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Converts reperation output to new data frame\
                                     containing specified information and saves it in gff3 format.')
    parser.add_argument("-i", "--inputTXT", action="store", dest="predictedORFs", required=True
                                          , help= "the input file. (created by reparation)")
    parser.add_argument("-c", "--condition", action="store", dest="condition", required=True
                                          , help= "the condition of the current file")
    parser.add_argument("-r", "--replicate", action="store", dest="replicate", required=True
                                           , help= "the replicate of the current file")
    parser.add_argument("-o", "--outputGFF", action="store", dest="outputGFF", required=True
                                           , help= "the output file name (gff3 format)")

    args = parser.parse_args()
    gff3df = to_gff3(args)

    gff3df.to_csv(args.outputGFF, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
