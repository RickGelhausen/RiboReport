rule retrieveGenome:
    input:
        "genome.fa"
    output:
        "genomes/genome.fa"
    threads: 1
    shell:
        "mkdir -p genomes; mv genome.fa genomes/"

rule retrieveAnnotation:
    input:
        "annotation.gtf"
    output:
        "annotation/annotation.gtf"
    threads: 1
    shell:
        "mkdir -p annotation; mv annotation.gtf annotation/"

# rule validateAnnotation:
#     input:
#         annotation="annotation.gtf"
#         genome="genome.fa"
#     output:
#         validAnnotation="annotation/annotation.gtf"
#         validGenome="genomes/genome.fa"
#     threads: 1
#     envs:
#         "../ribo_benchmark/envs/mergetools.yaml"
#     shell:
#         "mkdir -p annotation; mkdir -p genomes; ribo_benchmark/scripts/validateAnnotation -a {input.annotation} -g {input.genome} --annotationOutput {input.validAnnotation} --validGenome {input.validGenome}"
