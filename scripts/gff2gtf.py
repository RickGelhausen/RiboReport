#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections
import csv

def generate_annotation_dictionaries(annotation_path):
    """
    create dictionaries from annotation.
    gene_dict, cds_dict and rna_dict
    """

    annotation_df = pd.read_csv(annotation_path, sep="\t", comment="#", header=None)

    gene_dict = {}
    cds_dict = {}
    rna_dict = {}

    for row in annotation_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = int(getattr(row, "_3"))
        stop = int(getattr(row, "_4"))
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")

        attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            print(attribute_list)
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")

        if feature.lower() in ["gene","pseudogene"]:
            gene_id = ""
            if "id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("id")+1]

            gene_name = ""
            if "name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("name")+1]

            gene_biotype = ""
            if "gene_biotype" in attribute_list:
                gene_biotype = attribute_list[attribute_list.index("gene_biotype")+1]

            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

            gene_type = feature.lower()

            if gene_id == "":
                print("no valid ID for gene row.")
                print(row)

            counter = 1
            while gene_id in gene_dict:
                gene_id = "%s-%s" %(gene_id, counter)

                counter+=1
                
            gene_dict[gene_id] = (chromosome, start, stop, strand, gene_name, locus_tag, gene_biotype, gene_type)

        elif feature.lower() == "cds":
            parent = ""
            if "parent" in attribute_list:
                parent = attribute_list[attribute_list.index("parent")+1]

            cds_name = ""
            if "name" in attribute_list:
                cds_name = attribute_list[attribute_list.index("name")+1]

            if parent == "":
                print("no valid parent for cds row.")
                print(row)
            cds_dict[parent] = (cds_name)

        elif feature.lower() in ["rrna", "trna", "ncrna", "srna"]:
            parent = ""
            if "parent" in attribute_list:
                parent = attribute_list[attribute_list.index("parent")+1]

            if parent == "":
                print("no valid parent for rna row.")
                print(row)
            rna_dict[parent] = ()

    return gene_dict, cds_dict, rna_dict


def generate_gtf_file(in_gff, out_gtf):
    """
    read NCBI gff file, output ensembl gtf file
    """

    gene_dict, cds_dict, rna_dict = generate_annotation_dictionaries(in_gff)

    n_tuple = collections.namedtuple('Pandas', ["chrom", "source", "feature", "start", "stop", "score", "strand", "phase", "attributes"])
    counter = 0
    rows = []
    for idx, val in gene_dict.items():
        chrom, start, stop, strand, gene_name, locus_tag, gene_biotype, gene_type = val

        if locus_tag != "":
            gene_attributes = "gene_id \"%s\";" % locus_tag
        else:
            gene_attributes = "gene_id \"ltag%s\";" % counter

        if gene_name != "":
            gene_attributes += " gene_name \"%s\";" % gene_name
        else:
            gene_attributes += " gene_name \"gene%s\";" % counter

        if gene_biotype != "":
            gene_attributes += " gene_biotype \"%s\";" % gene_biotype
        else:
            if idx in cds_dict:
                gene_attributes += " gene_biotype \"protein_coding\";"
            elif idx in rna_dict:
                gene_attributes += " gene_biotype \"RNA\";"
            else:
                gene_attributes += " gene_biotype \"unknown\";"

        gene_row = [chrom, "converted", gene_type, start, stop, ".", strand, ".", gene_attributes]

        gene_attributes += " transcript_id \"transcript%s\";" % counter
        if gene_name != "":
            gene_attributes += " transcript_name \"%s-1\";" % gene_name
        else:
            gene_attributes += " transcript_name \"gene%s-1\";" % counter

        if gene_biotype != "":
            gene_attributes += " transcript_biotype \"%s\";" % gene_biotype
        else:
            if idx in cds_dict:
                gene_attributes += " transcript_biotype \"protein_coding\";"
            elif idx in rna_dict:
                gene_attributes += " transcript_biotype \"RNA\";"
            else:
                gene_attributes += " transcript_biotype \"unknown\";"

        transcript_row = [chrom, "converted", "transcript", start, stop, ".", strand, ".", gene_attributes]
        exon_row = [chrom, "converted", "exon", start, stop, ".", strand, ".", gene_attributes]

        rows.append(n_tuple(*gene_row))
        rows.append(n_tuple(*transcript_row))
        rows.append(n_tuple(*exon_row))

        if strand == "+":
            if idx in cds_dict:
                cds_row = [chrom, "converted", "CDS", start, stop-3, ".", strand, ".", gene_attributes]
            start_codon_row = [chrom, "converted", "start_codon", start, start+2, ".", strand, ".", gene_attributes]
            stop_codon_row = [chrom, "converted", "stop_codon", stop-2, stop, ".", strand, ".", gene_attributes]
        else:
            if idx in cds_dict:
                cds_row = [chrom, "converted", "CDS", start+3, stop, ".", strand, ".", gene_attributes]
            start_codon_row = [chrom, "converted", "start_codon", stop-2, stop, ".", strand, ".", gene_attributes]
            stop_codon_row = [chrom, "converted", "stop_codon", start, start+2, ".", strand, ".", gene_attributes]

        if idx in cds_dict:
            rows.append(n_tuple(*cds_row))
        rows.append(n_tuple(*start_codon_row))
        rows.append(n_tuple(*stop_codon_row))

        counter += 1

    return pd.DataFrame.from_records(rows, columns=["chrom", "source", "feature", "start", "stop", "score", "strand", "phase", "attributes"])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create gff2 in ensembl format from NCBI gff3.')
    parser.add_argument("-i", "--in_gff", action="store", dest="in_gff", required=True, help= "input gff file")
    parser.add_argument("-o", "--out_gtf", action="store", dest="out_gtf", required=True, help= "output gtf file")
    args = parser.parse_args()

    gtf_df = generate_gtf_file(args.in_gff, args.out_gtf)

    with open(args.out_gtf, "w") as f:
        gtf_df.to_csv(f, sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)




if __name__ == '__main__':
    main()
