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
    biotypeDict = dict()

    # biotype : chrom, start, stop, name, score, strand
    # biotype = ncRNA / cRNA
    for row in inputDF.itertuples(index=False, name='Pandas'):
        if getattr(row, "_2") == "gene":
            attributes = re.split('[;=]', getattr(row, "_8"))
            ID = ""
            if "ID" in attributes:
                ID = attributes[attributes.index("ID")+1]

            biotype = ""
            if "gene_biotype" in attributes:
                biotype = attributes[attributes.index("gene_biotype")+1]

            if biotype == "protein_coding":
                biotype = "coding"
            elif biotype in ["tRNA", "rRNA", "ncRNA"]:
                biotype = "noncoding"
            else:
                biotype = ""

            start = getattr(row, "_3")
            stop = getattr(row, "_4")
            strand = getattr(row, "_6")
            chromosome = getattr(row, "_0")

            if biotype in biotypeDict:
                biotypeDict[biotype].append((chromosome, start, stop, ID, "0", strand))
            else:
                biotypeDict[biotype] = []
    return biotypeDict


def create_bed(args):
    inputDF = pd.read_csv(args.inputGFF, sep='\t', header=None, comment="#")
    nTuple = collections.namedtuple('Pandas', ["chromosome","start","stop","name","score","strand"])

    # create a dictionary for common ids
    biotypeDict = create_dictionary(inputDF)

    number_noncoding = len(biotypeDict["noncoding"]) #+ 200
    # noncoding
    with open(args.outputBED+"_noncoding.bed", "w") as f:
        for chrom, start, stop, name, score, strand in biotypeDict["noncoding"]:
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start, stop, name, score, strand))

    # coding
    with open(args.outputBED+"_coding.bed", "w") as f:
        for chrom, start, stop, name, score, strand in biotypeDict["coding"][0:number_noncoding+1]:
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start, stop, name, score, strand))

    # all
    with open(args.outputBED+"_all.bed", "w") as f:
        for chrom, start, stop, name, score, strand in biotypeDict["coding"] + biotypeDict["noncoding"]:
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start, stop, name, score, strand, start, stop, "0,0,255", "1", stop-start+1, start))

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='removes duplicate intervals and finds\
                                                    the longest non-overlapping interval')
    parser.add_argument("-i", "--inputGFF", action="store", dest="inputGFF", required=True
                                          , help= "the input file (gff3 format).")
    parser.add_argument("-o", "--outputprefix", action="store", dest="outputBED", required=True
                                           , help= "the output file name (bed format)")
    args = parser.parse_args()
    if os.stat(args.inputGFF).st_size == 0:
       open(args.outputGFF, 'a').close()
    else:
       create_bed(args)

if __name__ == '__main__':
    main()
