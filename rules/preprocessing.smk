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

rule checkAnnotation:
    input:
        rules.retrieveAnnotation.output
    output:
        "annotation/annotation_processed.gff"
    threads: 1
    shell:
        "mkdir -p annotation; RiboReport/scripts/gtf2gff3.py -a {input} -o {output}"

rule createEnsemblAnnotation:
    input:
        rules.retrieveAnnotation.output
    output:
        "annotation/annotation_ensembl.gtf"
    threads: 1
    shell:
        "mkdir -p annotation; RiboReport/scripts/gff2gtf.py -i {input} -o {output}"