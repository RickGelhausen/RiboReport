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
        "annotation.gff"
    output:
        "annotation/annotation.gff"
    threads: 1
    shell:
        "mkdir -p annotation; mv annotation.gff annotation/"
