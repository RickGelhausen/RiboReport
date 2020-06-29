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

if TISHMODE == "RIBOONLY":
    rule all:
      input:
        expand("bam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("offsets/{method}-{condition}-{replicate}_p_offsets.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("reparation/{condition}-{replicate}/Predicted_ORFs.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"]),
        expand("coverage_deepribo/{condition}-{replicate}_cov_fwd.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("coverage_deepribo/{condition}-{replicate}_cov_rev.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("coverage_deepribo/{condition}-{replicate}_asite_fwd.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("coverage_deepribo/{condition}-{replicate}_asite_rev.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("deepribo/{condition}-{replicate}/data_list.csv", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("deepribo/{condition}-{replicate}/predictions.csv", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        "tracks/combined.gtf",
        "qc/multi/multiqc_report.html",
        expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/min/{method}-{condition}-{replicate}.min.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/min/{method}-{condition}-{replicate}.min.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/min/{method}-{condition}-{replicate}.min.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/min/{method}-{condition}-{replicate}.min.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("coverage/{method}-{condition}-{replicate}.bed", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        "tracks/potentialStopCodons.gff",
        "tracks/potentialStartCodons.gff",
        "tracks/potentialAlternativeStartCodons.gff",
        "tracks/potentialRibosomeBindingSite.gff",
        "auxiliary/final_annotation.xlsx",
        "auxiliary/final_annotation.gff",
        "auxiliary/final_annotation_complete.gff"
else:
    rule all:
      input:
        expand("bam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("offsets/{method}-{condition}-{replicate}_p_offsets.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("reparation/{condition}-{replicate}/Predicted_ORFs.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"]),
        expand("coverage_deepribo/{condition}-{replicate}_cov_fwd.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("coverage_deepribo/{condition}-{replicate}_cov_rev.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("coverage_deepribo/{condition}-{replicate}_asite_fwd.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("coverage_deepribo/{condition}-{replicate}_asite_rev.bedgraph", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("deepribo/{condition}-{replicate}/data_list.csv", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        expand("deepribo/{condition}-{replicate}/predictions.csv", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
        "tracks/combined.gtf",
        "qc/multi/multiqc_report.html",
        expand("transcripts/{condition}-{replicate}/transcripts.fa", zip, condition=samples.loc[samples["method"] == "RNA", "condition"], replicate=samples.loc[samples["method"] == "RNA", "replicate"]),
        expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/min/{method}-{condition}-{replicate}.min.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/min/{method}-{condition}-{replicate}.min.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/min/{method}-{condition}-{replicate}.min.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/min/{method}-{condition}-{replicate}.min.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("coverage/{method}-{condition}-{replicate}.bed", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        "tracks/potentialStopCodons.gff",
        "tracks/potentialStartCodons.gff",
        "tracks/potentialAlternativeStartCodons.gff",
        "tracks/potentialRibosomeBindingSite.gff",
        "auxiliary/final_annotation.xlsx",
        "auxiliary/final_annotation.gff",
        "auxiliary/final_annotation_complete.gff"


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

if TISHMODE != "RIBOONLY":
    include: "rules/irsom.smk"
    include: "rules/postprocessing_irsom.smk"
else:
    include: "rules/postprocessing.smk"

include: "rules/qc.smk"
include: "rules/visualization.smk"
include: "rules/readcounting.smk"
