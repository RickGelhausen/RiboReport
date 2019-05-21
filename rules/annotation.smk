rule ribotishAnnotation:
    input:
        annotation="annotation/annotation_togtf.gtf",
        sizes="genomes/sizes.genome"
    output:
        "annotation/annotation_processed.gtf"
    conda:
        "../envs/annotation.yaml"
    threads: 1
    shell:
        "mkdir -p annotation; ribo_benchmark/scripts/createRiboTISHannotation.py -a {input.annotation}  --genome_sizes {input.sizes} --annotation_output {output}"

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
