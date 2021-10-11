
rule predictSpectre:
   input:
       bam="bam/RIBO-{condition}-{replicate}.bam",
       bamindex="maplink/RIBO-{condition}-{replicate}.bam.bai",
       gtf="annotation/annotation_ensembl.gtf",
       fpkm="transcripts/{condition}-{replicate}/isoforms.fpkm_tracking",
   output:
       results="spectre/{condition}-{replicate}/result.txt",
       log="spectre/{condition}-{replicate}/spectre_log.txt"
   threads: 10
   singularity:
        "docker://gelhausr/spectre:latest"
   shell:
       """
       mkdir -p spectre;
       SPECtre.py --input {input.bam} --output {output.results} --log {output.log} --fpkm {input.fpkm} --gtf {input.gtf} --nt {threads}
       """

rule spectreGFF:
    input:
        "spectre/{condition}-{replicate}/result.txt"
    output:
        "spectre/{condition, [a-zA-Z]+}-{replicate,\d+}.spectre.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/spectreGFF.py -c {wildcards.condition}  -i {input} -o {output}"

rule concatSpectre:
    input:
        lambda wildcards: expand("spectre/{{condition}}-{replicate}.spectre.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.spectre.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/concatGFF.py {input} -o {output}"
