#!/usr/bin/env python
'''This script takes input files generated by
smorfer and creates a new data frame containing specified
information and writes it as gff3 format files.
'''

import pandas as pd
import argparse
import csv
import collections
import operator

# def create_nTuple(args, row):
#     nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])


#     chromosome = getattr(row, "_0")
#     start = str(getattr(row, "_1"))
#     stop = str(getattr(row, "_2"))
#     strand = getattr(row, "_5")

#     rpf = getattr(row, "_6")
#     no_clue = getattr(row, "_9")

#     source = "smorfer"
#     feature = "CDS"
#     score = rpf
#     phase = "."
#     attribute = "ID=" + chromosome + ":" + start + "-" + stop + ":" + strand \
#               + ";Name=" + chromosome + ":" + start + "-" + stop + ":" + strand \
#               + ";RPF=%s" % rpf + ";Condition=" + args.condition + ";Method=" + source

#     return nTuple(chromosome, source, feature, start, stop, score, strand, phase, attribute)


def create_stop_dict(args):
    in_df = pd.read_csv(args.in_file, sep='\t')
    in_df.columns = range(in_df.shape[1])

    stop_dict = {}
    for row in in_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        start = int(getattr(row, "_1"))
        stop = int(getattr(row, "_2"))
        strand = getattr(row, "_5")
        variance = getattr(row, "_6")
        rpf = getattr(row, "_8")
        if rpf == 0:
            continue
        if strand == "+":
            if (chromosome, stop, strand) in stop_dict:
                stop_dict[(chromosome, stop, strand)].append((variance, start, rpf))
            else:
                stop_dict[(chromosome, stop, strand)] = [(variance, start, rpf)]
        else:
            if (chromosome, start, strand) in stop_dict:
                stop_dict[(chromosome, start, strand)].append((variance, stop, rpf))
            else:
                stop_dict[(chromosome, start, strand)] = [(variance, stop, rpf)]

    return stop_dict

def to_gff3(args):
    stop_dict = create_stop_dict(args)
    nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])

    rows = []
    for key, value in stop_dict.items():
        chromosome, stop, strand = key
        variance, start, rpf = max(value, key=operator.itemgetter(0))

        source = "smorfer"
        feature = "CDS"
        score = rpf
        phase = "."

        if strand == "-":
            start, stop = stop, start

        start += 1

        attribute = "ID=" + chromosome + ":" + str(start) + "-" + str(stop) + ":" + strand \
                  + ";Name=" + chromosome + ":" + str(start) + "-" + str(stop) + ":" + strand \
                  + ";RPF=%s" % rpf + ";Condition=" + args.condition + ";Method=" + source

        rows.append(nTuple(chromosome, source, feature, start, stop, score, strand, phase, attribute))

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
