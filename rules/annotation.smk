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
        "annotation/annotation.gff"
    output:
        "auxiliary/enriched_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; RiboReport/scripts/enrich_annotation.py -a {input} -o {output}"

rule unambigousAnnotation:
    input:
        "auxiliary/enriched_annotation.gtf"
    output:
        "auxiliary/unambigous_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary;
        awk -F'\\t' '/^[^#]/ {{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\tID=uid%s;\\n", $1, $2, $3, $4, $5, $6, $7, $8, NR-1}}' {input} > {output}
        """
