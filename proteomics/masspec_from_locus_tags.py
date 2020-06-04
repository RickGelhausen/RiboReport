#!/usr/bin/env python
import os, sys
import pandas as pd
import argparse
import collections
import csv
import re

def create_dictionary(args):
    gff_df = pd.read_csv(args.input_gff, sep="\t", comment="#", header=None)

    gene_dict = {}
    for row in gff_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")

        attribute_list = [x for x in re.split('[;=]', attributes)]
        if feature.lower() == "gene":
            gene_id = attribute_list[attribute_list.index("ID") + 1]
            if "old_locus_tag" in attributes:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag") + 1]
                gene_dict[gene_id] = old_locus_tag

    return gene_dict

def get_massspec_data(args, gene_dict, locus_tags):
    gff_df = pd.read_csv(args.input_gff, sep="\t", comment="#", header=None)

    rows = []
    old = []
    for row in gff_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")

        attribute_list = [x for x in re.split('[;=]', attributes)]
        if feature.lower() == "pseudogene":
            if "old_locus_tag" not in attributes:
                continue
            old_locus_tag = attribute_list[attribute_list.index("old_locus_tag") + 1]
            if old_locus_tag in locus_tags:
                old.append(old_locus_tag)
                rows.append(row)
            else:
                continue
        elif feature.lower() != "gene":
            if "Parent" in attributes:
                parent = attribute_list[attribute_list.index("Parent") + 1]
                if parent not in gene_dict:
                    continue
                if gene_dict[parent] in locus_tags:
                    old.append(gene_dict[parent])
                    rows.append(row)
                else:
                    continue
            else:
                continue
        else:
            continue

    print(len([ x for x in locus_tags if x not in old]))

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

def read_tsv(args):
    return list(pd.read_csv(args.massspec, sep="\t", comment="#", header=None)[0])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Transfer the labeled data, to the complete gff file')
    # required
    parser.add_argument("-a", "--annotation", action="store", dest="input_gff", required=True, help="the annotation to become labeled")
    parser.add_argument("-m", "--massspec", action="store", dest="massspec", required=True, help="list of massspec locus tags")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help="the output file")
    args = parser.parse_args()

    out_df = get_massspec_data(args, create_dictionary(args), read_tsv(args))
    out_df.to_csv(args.output, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
