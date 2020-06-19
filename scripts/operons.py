#!/usr/bin/env python
# time ./operons.py --in_gff_filepath ecoli/annotation.gff
import argparse
import os
import csv as csv
import interlap
from operator import attrgetter


class Feature:
    def __init__(self, gfftype, seqid, start, end, strand, extended_start, extended_end, attributes):
        self.gfftype = gfftype
        self.seqid = seqid
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.extended_start = int(extended_start)
        self.extended_end = int(extended_end)
        self.attributes = attributes
    def __str__(self):
        return self.gfftype + "_" + self.seqid + "_" + str(self.start) + "_" + str(self.end) + "_" + self.strand + "_" + str(self.extended_start) + "_" + str(self.extended_end) + "_" + self.attributes + "\n" 
    def gff(self):
        return self.seqid + "\t.\t" + self.gfftype + "\t" + str(self.start) + "\t" + str(self.end) + "\t.\t" + self.strand + "\t.\t" + "ID=operon" + str(self.start) + "_" + str(self.end) + ";" + self.attributes + "\n"

def get_overlapping_cds(cds_key,cdss,inter,operon_cds,detected_keys):
    current_cds = cdss[cds_key]
    cds_interval = (int(current_cds.extended_start), int(current_cds.extended_end))
    overlapping_intervals = list(inter.find(cds_interval))
    if not len(overlapping_intervals) == 1:
        operon_cds.append(current_cds)
        for cds_interval2 in overlapping_intervals:
            #intersectiv = get_overlap_bounderies(cds_interval, cds_interval2)
            operon_cds_key = str(cds_interval2[0]) + "-" + str(cds_interval2[1])
            if not operon_cds_key in detected_keys:
                if operon_cds_key in cdss.keys():
                    detected_keys[operon_cds_key]=""
                    #operon_cds.append(operon_cds_key)
                    #delete key from dict
                    #del cdss[operon_cds_key]
                    get_overlapping_cds(operon_cds_key,cdss,inter,operon_cds,detected_keys)
                    #continue here
    return operon_cds

def operons(seqid, cdss, strand):
    inter = interlap.InterLap()
    intervals = []
    operons = []
    for cds_key in cdss:
        current_cds = cdss[cds_key]
        if current_cds.seqid == seqid:
            intervals.append((int(current_cds.extended_start), int(current_cds.extended_end)))
    inter.update(intervals)
    detected_keys = {}
    for cds2_key in cdss:
        if not cdss[cds2_key].seqid == seqid:
            continue
        if not cds2_key in detected_keys:
            current_cds2 = cdss[cds2_key]
            operon_cds = []
            overlapping_cds = get_overlapping_cds(cds2_key,cdss,inter, operon_cds, detected_keys)
            attributes= "N_cds=" + str(len(overlapping_cds)) + ";"
            if len(overlapping_cds) == 0:
                continue
            else:
                #operon_start = str(min(cds.start for cds in overlapping_cds))
                cds_starts = []
                for cds3 in overlapping_cds:
                #    print("Starts")
                #    print(cds3.start)
                    cds_starts.append(cds3.start)
                #print("Operon_start")
                operon_start = min(cds_starts)
                #print(operon_start)
                #operon_end = str(max(cds.end for cds in overlapping_cds))
                operon_end = str(max(overlapping_cds,key=attrgetter('end')))
                cds_end = []
                for cds3 in overlapping_cds:
                #    print("Ends")
                #    print(cds3.end)
                    cds_end.append(cds3.end)
                #print("Operon_end")
                operon_end = max(cds_end)
                #print(operon_end)

                new_operon = Feature("operon", seqid, operon_start, operon_end, strand, str(0), str(0), attributes)
                #debug
                #print("Operon") 
                #print(new_operon)
                #print("CDS_list")
                #overlapping_cds_string=''.join([str(x) for x in overlapping_cds])
                #print(overlapping_cds_string)
                operons.append(new_operon)
    return operons

def get_cds(input_gff_filepath):
    cdsminus = {}
    cdsplus = {}
    seqid_set = set()
    with open(input_gff_filepath, newline='\n') as csvfile:
        gffreader = csv.reader(csvfile, delimiter='\t')
        for entry in gffreader:
            if (not entry[0].startswith('#')):
                gfftype = entry[2]
                if gfftype == "CDS":
                    seqid = entry[0]
                    seqid_set.add(seqid)
                    source = entry[1]
                    gfftype = entry[2]
                    start = entry[3]
                    end = entry[4]
                    score = entry[5]
                    strand = entry[6]
                    phase = entry[7]
                    attributes = entry[8]
                    if strand == "+":
                        extended_start = str(int(start) - 100)
                        extended_end = str(int(end) + 100)
                        key = str(extended_start) + "-" + str(extended_end)
                        cdsplus[key]=Feature(gfftype, seqid, start, end, strand, extended_start, extended_end,"")
                    else:
                        extended_start = str(int(start) - 100)
                        extended_end = str(int(end) + 100)
                        key = str(extended_start) + "-" + str(extended_end)
                        cdsminus[key]=Feature(gfftype, seqid, start, end, strand, extended_start, extended_end,"")
    seqids = list(seqid_set)
    return (seqids, cdsminus, cdsplus)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Operons.')
    parser.add_argument("--in_gff_filepath", help='Input gff path', required=True)
    args = parser.parse_args()
    (seqids,cdsminus,cdsplus)=get_cds(args.in_gff_filepath)
    #print(cds)
    strands = ["+","-"]
    operons_list = [] 
    for seqid in seqids:
        minus_operons = operons(seqid, cdsminus, "-")
        operons_list.extend(minus_operons)
        plus_operons = operons(seqid, cdsplus, "+")
        operons_list.extend(plus_operons)
        operons_string=''.join([str(x.gff()) for x in operons_list])
        print(operons_string)

        #print(operons_list)
if __name__ == '__main__':
    main()
