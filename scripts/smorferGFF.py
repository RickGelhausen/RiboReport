#!/usr/bin/env python
'''This script takes input files generated by
smorfer and creates a new data frame containing specified
information and writes it as gff3 format files.
'''

import pandas as pd
import argparse
import csv
import collections


def create_nTuple(args, row):
    nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])


    chromosome = getattr(row, "_0")
    start = str(getattr(row, "_1"))
    stop = str(getattr(row, "_2"))
    strand = getattr(row, "_5")

    rpf = getattr(row, "_6")
    no_clue = getattr(row, "_9")

    source = "smorfer"
    feature = "CDS"
    score = rpf
    phase = "."
    attribute = "ID=" + chromosome + ":" + start + "-" + stop + ":" + strand \
              + ";Name=" + chromosome + ":" + start + "-" + stop + ":" + strand \
              + ";RPF=%s" % rpf + ";Condition=" + args.condition + ";Method=" + source

    return nTuple(chromosome, source, feature, start, stop, score, strand, phase, attribute)


def to_gff3(args):
    in_df = pd.read_csv(args.in_file, sep='\t', header=None)

    rows = []
    for row in in_df.itertuples(index=False, name='Pandas'):
        rows.append(create_nTuple(args,row))
    rows = [row for row in rows if row is not None]

    res_df = pd.DataFrame.from_records(rows, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"])

    return res_df


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Converts smorfer output to new data frame\
                                     containing specified information and saves it in gff3 format.')
    parser.add_argument("-i", action="store", dest="in_file", required=True
                                          , help= "the unfiltered tsv file (created by price)")
    parser.add_argument("-c", action="store", dest="condition", required=True
                                          , help= "the condition of the current file")
    parser.add_argument("-o", action="store", dest="out_file", required=True
                                           , help= "the output file name (gff3 format)")
    args = parser.parse_args()
    gff3_df = to_gff3(args)

    gff3_df.to_csv(args.out_file, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
