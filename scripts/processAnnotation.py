import csv
import pandas as pd
import argparse

def process_annotation(args):
    # read input annotation
    annDF = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)

    # check if gff2/gtf or gff3 format
    if annDF[8].str.contains("ID="):
        # gff3
        annDF.to_csv(args.output, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    else:
        # only accept rows that contain gene_id transcript_id combination
        annDF = annDF[annDF[8].str.contains("gene_id") & annDF[8].str.contains("transcript_id")]

        # write processed annotation to output
        annDF.to_csv(args.output, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='process the annotation for use with plastid')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "the input annotation file.")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "the processed annotation file.")
    args = parser.parse_args()

    process_annotation(args)


if __name__ == '__main__':
    main()
