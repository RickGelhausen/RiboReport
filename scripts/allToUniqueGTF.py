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
        
        if "Condition" in attributes and "Method" in attributes:
            condition = attributes[attributes.index("Condition")+1]
            method = attributes[attributes.index("Method")+1]

        # (id, tool, condition) = (start, stop, strand, score)
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")

        if getattr(row, "_1") == "ribotish":
            score = float(attributes[attributes.index("Ribo_pvalue")+1])
        elif getattr(row, "_1") == "reparation":
            score = 1 - float(attributes[attributes.index("Prob")+1])
        elif getattr(row, "_1") == "deepribo":
            score = 0.0001 # take all
        elif getattr(row, "_1") == "irsom":
            score = 1 - float(attributes[attributes.index("Prob")+1])

        idx = (getattr(row, "_0"), method, condition, start, stop, strand)
        if idx in geneDict:
            if geneDict[idx] > score:
                geneDict[idx] = score
        else:
            geneDict[idx] = score
    return geneDict

def create_gtf(args):
    inputDF = pd.read_csv(args.inputGFF, sep='\t', header=None)
    nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])

    # create a dictionary for common ids
    geneDict = create_dictionary(inputDF)

    # run over all entries in the dictionary and combine overlapping ones
    rows = []
    for key in geneDict.keys():
        score = geneDict[key]
        accession = key[0]
        method = key[1]
        condition = key[2]
        start = key[3]
        stop = key[4]
        strand = key[5]

        # new content
        seqName = "1"
        source = method
        type = "CDS"
        phase = "."
        id = "%s:%s-%s:%s" % (accession, start, stop, strand)

        attribute = "gene_id \"%s\"; method \"%s\"; condition \"%s\"" % (id, method, condition)
        rows.append(nTuple(seqName, source, type, start, stop, score, strand, phase, attribute))

    return pd.DataFrame.from_records(rows, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='removes duplicate intervals and finds\
                                                    the longest non-overlapping interval')
    parser.add_argument("-i", "--inputGFF", action="store", dest="inputGFF", required=True
                                          , help= "the input file (gff3 format).")
    parser.add_argument("-o", "--outputGTF", action="store", dest="outputGTF", required=True
                                           , help= "the output file name (gtf format)")
    args = parser.parse_args()
    if os.stat(args.inputGFF).st_size == 0:
       open(args.outputGFF, 'a').close()
    else:
       newDF = create_gtf(args)
       newDF.to_csv(args.outputGTF, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
