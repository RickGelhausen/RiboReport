rule ribotishAnnotation:
    input:
        annotation="qc/featurecount/annotation.gtf",
        sizes="genomes/sizes.genome"
    output:
        "annotation/annotation_processed.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p ribotish; RiboReport/scripts/createRiboTISHannotation.py -a {input.annotation} --genome_sizes {input.sizes} --annotation_output {output}"

rule gff3ToGenePred:
    input:
        "annotation/annotation.gff"
    output:
        "annotation/annotation.genepred"
    conda:
        "../envs/annotation.yaml"
    threads: 1
    shell:
        "mkdir -p ribotish; gff3ToGenePred {input} {output}"

rule genePredToBed:
    input:
        "annotation/annotation.genepred"
    output:
        "annotation/annotation.bed"
    conda:
        "../envs/annotation.yaml"
    threads: 1
    shell:
        "mkdir -p ribotish; genePredToBed {input} {output}"

rule featurecountAnnotation:
    input:
        annotation={rules.retrieveAnnotation.output},
    output:
        "qc/all/annotationall.gtf",
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p qc/all; RiboReport/scripts/annotation_featurecount.py -a {input.annotation} -o {output};"

rule enrichAnnotation:
    input:
        annotation=rules.checkAnnotation.output
    output:
        "auxiliary/enriched_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; RiboReport/scripts/enrich_annotation.py -a {input.annotation} -o {output}"

rule unambigousAnnotation:
    input:
        "auxiliary/enriched_annotation.gff"
    output:
        "auxiliary/unambigous_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary;
        awk -F'\\t' '/^[^#]/ {{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\tID=uid%s;\\n", $1, $2, $3, $4, $5, $6, $7, $8, NR-1}}' {input} > {output}
        """

rule generateTranscripts:
    input:
        bam="maplink/RNA-{condition}-{replicate}.bam",
        bamindex="maplink/RNA-{condition}-{replicate}.bam.bai",
        annotation=rules.checkAnnotation.output
    output:
        transcripts="transcripts/{condition}-{replicate}/transcripts.gtf",
        isoforms="transcripts/{condition}-{replicate}/isoforms.fpkm_tracking",
        genes="transcripts/{condition}-{replicate}/genes.fpkm_tracking",
    conda:
        "../envs/cufflinks.yaml"
    threads: 20
    shell:
        "mkdir -p transcripts; cufflinks {input.bam} -p {threads} -o ./transcripts/{wildcards.condition}-{wildcards.replicate}/ -g {input.annotation} --library-type fr-firststrand"


