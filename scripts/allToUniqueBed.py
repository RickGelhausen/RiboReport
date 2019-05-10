#!/usr/bin/env python
'''This script takes input gff3 files and handles
overlapping intervals, by removing duplicates and
finding the longest non-overlapping interval.
'''
from operator import itemgetter
import pandas as pd
import re
import argparse
import numpy as np
import os
import csv
import collections

# helper function to create a dictionary {geneID : namedtuple}
def create_dictionary(inputDF):
    geneDict = dict()
    for row in inputDF.itertuples(index=False, name='Pandas'):
        attributes = re.split('[;=]', getattr(row, "_8"))
        # if "ID" in attributes:
        #      geneID = attributes[attributes.index("ID")+1]
        if "Condition" in attributes and "Method" in attributes:
            condition = attributes[attributes.index("Condition")+1]
            method = attributes[attributes.index("Method")+1]

        # (id, tool, condition) = (start, stop, strand, score)
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")

        if getattr(row, "_1") == "ribotish":
            score = attributes[attributes.index("Ribo_pvalue")+1]
        elif getattr(row, "_1") == "reparation":
            score = attributes[attributes.index("Prob")+1]
        elif getattr(row, "_1") == "deepribo":
            score = attributes[attributes.index("pvalue")+1]


        idx = (getattr(row, "_0"), method, condition, start, stop, strand)
        if idx in geneDict:
            geneDict[idx].append(score)
        else:
            geneDict[geneID] = [score]
    return geneDict

def create_bed(args):
    inputDF = pd.read_csv(args.inputGFF, sep='\t', header=None)

    # create a dictionary for common ids
    geneDict = create_dictionary(inputDF)

    # run over all entries in the dictionary and combine overlapping ones
    rows = []
    for key in geneDict.keys():
        scores = geneDict[key]
        rows.append(nTuple(key[0], key[1], key[2], key[3], keys[4], keys[5], ",".join(scores)))

    return pd.DataFrame.from_records(rows, columns=["id", "method", "condition", "start", "stop", "strand", "scores"])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='removes duplicate intervals and finds\
                                                    the longest non-overlapping interval')
    parser.add_argument("-i", "--inputGFF", action="store", dest="inputGFF", required=True
                                          , help= "the input file (gff3 format).")
    parser.add_argument("-o", "--outputBED", action="store", dest="outputBED", required=True
                                           , help= "the output file name (bed format)")
    args = parser.parse_args()
    if os.stat(args.inputGFF).st_size == 0:
       open(args.outputGFF, 'a').close()
    else:
       newDF = create_bed(args)
       newDF.to_csv(args.outputBED, sep="\t", header=True, index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
