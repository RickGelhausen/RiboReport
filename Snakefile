import os
import re
import pandas as pd
import itertools as iter
from snakemake.utils import validate, min_version
min_version("5.4.2")

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
    expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"])

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
include: "rules/auxillary.smk"
#indexing
include: "rules/indexing.smk"
if TISHMODE == "RIBOONLY":
   #ribotish
   include: "rules/ribotish.smk"
else:
   #ribotish
   include: "rules/ribotishall.smk"
#reparation
include: "rules/reparation.smk"
