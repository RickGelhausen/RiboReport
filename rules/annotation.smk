# rule processAnnotation:
#     input:
#         annotation=rules.retrieveAnnotation.output
#     output:
#         "annotation/processed-annotation.gtf"
#     conda:
#         "../envs/pytools.yaml"
#     threads: 1
#     shell:
#         "mkdir -p annotation; python3 ribo_benchmark/scripts/processAnnotation.py -a {input.annotation} -o {output}"

rule ribotishAnnotation:
    input:
        annotation="auxiliary/featurecount/annotation.gtf",
        sizes="genomes/sizes.genome"
    output:
        "annotation/annotation_processed.gtf"
    conda:
        "../envs/annotation.yaml"
    threads: 1
    shell:
        "mkdir -p annotation; ribo_benchmark/scripts/createRiboTISHannotation.py -a {input.annotation}  --genome_sizes {input.sizes} --annotation_output {output}"

#rule ribocodeAnnotation:
#    input:
#        annotation="auxiliary/featurecount/annotation.gtf",
#        sizes="genomes/sizes.genome"
#    output:
#        "ribocode/annotation.gtf"
#    conda:
#        "../envs/annotation.yaml"
#    threads: 1
#    shell:
#        "mkdir -p annotation; ribo_benchmark/scripts/createRiboCodeAnnotation.py -a {input.annotation}  --genome_sizes {input.sizes} --annotation_output {output}"

rule ribocodeAnnotation:
    input:
        "annotation/annotation_ensembl.gtf"
    output:
        "ribocode/annotation.gtf"
    conda:
        "../envs/ribocode.yaml"
    threads: 1
    shell:
        "mkdir -p annotation; GTFupdate {input} > {output}"


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
