import os
import re
import pandas as pd
import itertools as iter
from snakemake.utils import validate, min_version
min_version("5.4.5")

ADAPTERS=config["adapter"]
INDEXPATH=config["genomeindexpath"]
CODONS=config["alternativestartcodons"]
TISHMODE=config["tishmode"]

onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_csv(config["samples"], dtype=str, sep="\t").set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])
validate(samples, schema="schemas/samples.schema.yaml")

rule all:
  input:
    expand("bam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
    expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
    expand("offsets/{method}-{condition}-{replicate}_p_offsets.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
    expand("reparation/{condition}-{replicate}/Predicted_ORFs.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
    expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"]),
    expand("coverage/{condition}-{replicate}_cov_fwd.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
    expand("coverage/{condition}-{replicate}_cov_rev.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
    expand("coverage/{condition}-{replicate}_asite_fwd.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
    expand("coverage/{condition}-{replicate}_asite_rev.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
    expand("deepribo/{condition}-{replicate}/data_list.csv", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
    expand("deepribo/{condition}-{replicate}/predictions.csv", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
    "tracks/combined.gtf",
    "qc/multi/multiqc_report.html",


    #expand("irsom/{condition}-{replicate}/test.txt", zip, condition=samples.loc[samples["method"] == "RNA", "condition"], replicate=samples.loc[samples["method"] == "RNA", "replicate"]),


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
#coverage
include: "rules/coverage.smk"
#auxillary
include: "rules/auxiliary.smk"
#indexing
include: "rules/indexing.smk"
#ribotish
include: "rules/annotation.smk"
if TISHMODE == "RIBOONLY":
   include: "rules/ribotish.smk"
else:
   include: "rules/ribotishall.smk"
#reparation
include: "rules/reparation.smk"
#deepribo
include: "rules/deepribo.smk"
#rnacode
#include: "rules/rnacode.smk"
include: "rules/postprocessing.smk"
include: "rules/irsom.smk"
include: "rules/qc.smk"
