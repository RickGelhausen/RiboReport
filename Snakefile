import os
import re
import pandas as pd
import itertools as iter
from snakemake.utils import validate, min_version
min_version("5.4.5")

ADAPTERS=config["adapter"]
CODONS=config["alternativestartcodons"]

onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_csv(config["samples"], dtype=str, sep="\t").set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])
validate(samples, schema="schemas/samples.schema.yaml")

hasRIBO=True
if "RIBO" not in samples["method"].unique():
    hasRIBO=False
    print("No Ribo-seq libraries were detected. No prediction tools for this setup are currently implemented. If you have pure Ribo-seq libraries, please use the method tag RIBO. Continuing...")

def get_wigfiles(wildcards):
  method=samples["method"]
  condition=samples["condition"]
  replicate=samples["replicate"]
  wilds = zip(method, condition, replicate)

  bigwigs = [["global", "centered", "fiveprime", "threeprime"], ["raw", "mil", "min"], ["forward", "reverse"], list(wilds)]
  bigwigs = list(iter.product(*bigwigs))

  wigfiles = []
  for bw in bigwigs:
      wigfiles.append("%stracks/%s/%s-%s-%s.%s.%s.%s.bw" %(bw[0], bw[1], bw[3][0], bw[3][1], bw[3][2], bw[1], bw[2], bw[0]))

  return wigfiles

if hasRIBO:
    rule all:
      input:
        expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("reparation/{condition}-{replicate}/Predicted_ORFs.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"]),
        expand("deepribo/{condition}-{replicate}/predictions.csv", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        "qc/multi/multiqc_report.html",
        get_wigfiles,
        "tracks/potentialStopCodons.gff",
        "tracks/potentialStartCodons.gff",
        "tracks/potentialAlternativeStartCodons.gff",
        "tracks/potentialRibosomeBindingSite.gff",
        "auxiliary/final_annotation.xlsx",
        "auxiliary/final_annotation.gff",
        "auxiliary/final_annotation_complete.gff",
        "tracks/predictions.gtf"

else:
    print("No Ribo libraries given")



onsuccess:
    print("Done, no error")

#Preprocessing
include: "rules/preprocessing.smk"
#Adaper removal and quality control
include: "rules/trimming.smk"
#removal of reads mapping to ribosomal rna genes
include: "rules/rrnafiltering.smk"
#mapping
include: "rules/mapping.smk"
#maplink
include: "rules/maplink.smk"
include: "rules/maplinktis.smk"
#auxillary
include: "rules/auxiliary.smk"
#indexing
include: "rules/indexing.smk"
#ribotish
include: "rules/annotation.smk"
include: "rules/ribotish.smk"

#reparation
include: "rules/reparation.smk"
#deepribo
include: "rules/deepribo.smk"
#rnacode
#include: "rules/rnacode.smk"

include: "rules/irsom.smk"
include: "rules/postprocessing.smk"
include: "rules/spectre.smk"
include: "rules/qc.smk"
include: "rules/visualization.smk"
include: "rules/readcounting.smk"
