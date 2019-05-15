rule generateTranscripts:
    input:
        bam="maplink/RNA-{condition}-{replicate}.bam",
        bamindex="maplink/RNA-{condition}-{replicate}.bam.bai"
    output:
        "transcripts/{condition}-{replicate}/transcripts.gtf"
    conda:
        "../envs/cufflinks.yaml"
    threads: 20
    shell:
        "mkdir -p transcripts; cufflinks {input.bam} -p {threads} -o ./transcripts/{wildcards.condition}-{wildcards.replicate}/"

rule createCodingBed:
    input:
        "annotation/annotation.gtf"
    output:
        coding="cpat/annotation_coding.bed",
        noncoding="cpat/annotation_noncoding.bed",
        all="cpat/annotation_all.bed"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p cpat; ribo_benchmark/scripts/retrieveCodingInformation.py -i {input} -o cpat/annotation"

rule bedToFasta:
    input:
        genome="genomes/genome.fa",
        coding="cpat/annotation_coding.bed",
        noncoding="cpat/annotation_noncoding.bed"
    output:
        coding="cpat/annotation_coding.fa",
        noncoding="cpat/annotation_noncoding.fa"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        """
        mkdir -p cpat;
        bedtools getfasta -fi {input.genome} -bed {input.coding} -fo {output.coding}
        bedtools getfasta -fi {input.genome} -bed {input.noncoding} -fo {output.noncoding}
        """

rule makeHexamerTab:
    input:
        coding="cpat/annotation_coding.fa",
        noncoding="cpat/annotation_noncoding.fa"
    output:
        "cpat/hexamer.tsv"
    threads: 1
    shell:
        """
        source activate cpat_fix; 
        mkdir -p cpat; 
        make_hexamer_tab.py -c {input.coding} -n {input.noncoding} > {output};
        """

rule makeLogitModel:
    input:
        hexamer="cpat/hexamer.tsv",
        coding="cpat/annotation_coding.fa",
        noncoding="cpat/annotation_noncoding.fa"
    output:
          xls="cpat/model.feature.xls",
          model="cpat/model.make_logitModel.r"
    threads: 1
    shell:
        """
        source activate cpat_fix;
        make_logitModel.py -x {input.hexamer} -c {input.coding} -n {input.noncoding} -o cpat/model;
        """

rule logit:
    input:
        xls="cpat/model.feature.xls",
        model="cpat/model.make_logitModel.r"
    output:
        "cpat/model.logit.RData"
    threads: 1
    conda:
        "../envs/rtools.yaml"
    shell:
        "mkdir -p cpat; Rscript {input.model}"


rule transcriptsBED:
    input:
        "transcripts/{condition}-{replicate}/transcripts.gtf"
    output:
        "transcripts/{¢ondition}-{replicate}/transcripts.bed"
    threads: 1
    conda:
        "../envs/mergetools.yaml"
    shell:
        "mkdir -p transcripts; ribo_benchmark/scripts/transcriptsToBed.py -i {input} -o {output}"

rule cpat:
    input:
        genome="genomes/genome.fa",
        ref="transcripts/{condition}-{replicate}/transcripts.bed",
        hexamer="cpat/hexamer.tsv",
        logit="cpat/model.logit.RData"
    output:
        script="cpat/{condition}-{replicate}.r",
        dat="cpat/{condition}-{replicate}.dat"
    threads: 1
    shell:
        """
        source activate cpat_fix;
        cpat.py -r {input.genome} -g {input.annotation} -x {input.hexamer} -d {input.logit} -o cpat/{wildcards.condition}-{wildcards.replicate}
        """

rule cpatScript:
    input:
        xls="cpat/{condition}-{replicate}.dat",
        model="cpat/{condition-{replicate}.r"
    output:
        "cpat/{condition}-{replicate}.tsv"
    threads: 1
    conda:
        "../envs/rtools.yaml"
    shell:
        "mkdir -p cpat; Rscript {input.model}; mv cpat/{wildcards.condition}-{wildcards.replicate} cpat/{wildcards.condition}-{wildcards.replicate}.tsv"

