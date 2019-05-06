rule prepareTranscripts:
    input:
        annotation="ribocode/annotation.gtf",
        genome="genomes/genome.fa"
    output:
        "ribocode/annotation/transcripts_cds.txt"
    conda:
        "../envs/ribocode.yaml"
    threads: 1
    shell:
        "mkdir -p ribocode; prepare_transcripts -g {input.annotation} -f {input.genome} -o ribocode/annotation/"

rule selectingRPFlength:
    input:
        annotation="ribocode/annotation/transcripts_cds.txt",
        bam="maplink/RIBO-{condition}-{replicate}.bam"
    output:
        "ribocode/{condition}-{replicate}/_pre_config.txt"
    conda:
        "../envs/ribocode.yaml"
    threads: 1
    shell:
        "mkdir -p ribocode; metaplots -a ribocode/annotation/ -r {input.bam} -o ribocode/{wildcards.condition}-{wildcards.replicate}/"

rule ribocode:
    input:
        annotation="ribocode/annotation/transcripts_cds.txt",
        config="ribocode/{condition}-{replicate}/_pre_config.txt",
    output:
       "ribocode/{condition}-{replicate}/ribocode_orfs.txt"
    conda:
        "../envs/ribocode.yaml"
    threads: 1 
    shell:
        "mkdir -p ribocode; RiboCode -a ribocode/annotation/ -c {input.config} -l no -g -o {output}"
