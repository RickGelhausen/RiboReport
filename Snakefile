import os
import re
import pandas as pd
import itertools as iter
from snakemake.utils import validate, min_version
min_version("5.4.2")

ADAPTERS=config["adapter"]
INDEXPATH=config["genomeindexpath"]
CODONS=config["alternativestartcodons"]

onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_csv(config["samples"], dtype=str, sep="\t").set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])
validate(samples, schema="schemas/samples.schema.yaml")


rule all:
  input:
    expand("bam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])

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
