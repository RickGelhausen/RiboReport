#!/usr/bin/env python
'''This script takes input files generated by
deepribo and creates a new data frame containing specified
information and writes it as gff3 format files.
'''

import pandas as pd
import argparse
import csv
import collections


# little helper function to create named tuple without having to always state every argument
def createNTuple(args, row):
    nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])
    # txt file content
    filename = getattr(row, "filename")
    filename_counts = getattr(row, "filename_counts")
    label = bool(getattr(row, "label"))
    in_gene = getattr(row, "in_gene")
    strand = str(getattr(row, "strand"))
    coverage = str(getattr(row, "coverage"))
    coverage_elo = str(getattr(row, "coverage_elo"))
    rpk = str(getattr(row, "rpk"))
    rpk_elo = str(getattr(row, "rpk_elo"))
    start_site = str(getattr(row, "start_site"))
    start_codon = getattr(row, "start_codon")
    stop_site = str(getattr(row, "stop_site"))
    stop_codon = str(getattr(row, "stop_codon"))
    locus = str(getattr(row, "locus"))
    prot_seq = str(getattr(row, "prot_seq"))
    nuc_seq = str(getattr(row, "nuc_seq"))
    pred = str(getattr(row, "pred"))
    pred_rank = str(getattr(row, "pred_rank"))
    SS = str(getattr(row, "SS"))
    dist = str(getattr(row, "dist"))
    SS_pred_rank = getattr(row, "SS_pred_rank")

    # new content
    chromosome, rest = locus.split(":")
    start, stop = rest.split("-")

    if SS_pred_rank == 999999:
        return

    seqName = chromosome
    source = "deepribo"
    feature = "CDS"
    score = "."
    phase = "."
    attribute = "ID=" + chromosome + ":" + start + "-" + stop + ":" + strand \
              + ";Name=" + chromosome + ":" + start + "-" + stop + ":" + strand \
              + ";SS_pred_rank=" + str(SS_pred_rank) + ";Condition=" + args.condition + ";Method=deepribo"

    return nTuple(seqName, source, feature, start, stop, score, strand, phase, attribute)


def to_gff3(args):
    inputDF = pd.read_csv(args.predictedORFs, sep=',')

    # extract information from each row and build new dataframe in gff format
    rows = []
    for row in inputDF.itertuples(index=True, name='Pandas'):
        rows.append(createNTuple(args, row))
    rows = [row for row in rows if row is not None]
    return pd.DataFrame.from_records(rows, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"])


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Converts deepribo output to new data frame\
                                     containing specified information and saves it in gff3 format.')
    parser.add_argument("-i", "--inputCSV", action="store", dest="predictedORFs", required=True
                                          , help= "the input file. (created by reparation)")
    parser.add_argument("-c", "--condition", action="store", dest="condition", required=True
                                          , help= "the condition of the current file")
    parser.add_argument("-o", "--outputGFF", action="store", dest="outputGFF", required=True
                                           , help= "the output file name (gff3 format)")

    args = parser.parse_args()
    gff3df = to_gff3(args)

    gff3df.to_csv(args.outputGFF, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
