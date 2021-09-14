#!/usr/bin/env python
'''This script takes input gtf files and handles
and converts it to a bed12 file
'''

import pandas as pd
import re
import argparse

import os


def create_bed(args):
    inputDF = pd.read_csv(args.inputGTF, sep='\t', header=None)
    with open(args.outputBED, "w") as out:
        for row in inputDF.itertuples(index=False, name='Pandas'):
            if getattr(row, "_2").lower() == "transcript":
                attributes = re.split('[; ]', getattr(row, "_8"))

                start = getattr(row, "_3")
                stop = getattr(row, "_4")
                strand = getattr(row, "_6")
                chromosome = getattr(row, "_0")
                name = "%s:%s-%s:%s" % (chromosome, start, stop, strand)

                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chromosome, start, stop, name, "0", strand, start, stop, "0,0,255", "1", stop-start+1, start))

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='removes duplicate intervals and finds\
                                                    the longest non-overlapping interval')
    parser.add_argument("-i", "--inputGTF", action="store", dest="inputGTF", required=True
                                          , help= "the input file (gtf format).")
    parser.add_argument("-o", "--outputBED", action="store", dest="outputBED", required=True
                                           , help= "the output file name (bed12 format)")
    args = parser.parse_args()
    if os.stat(args.inputGTF).st_size == 0:
       open(args.outputBED, 'a').close()
    else:
       create_bed(args)

if __name__ == '__main__':
    main()
