rule generateTranscriptsSpectre:
    input:
        bam="maplink/RNA-{condition}-{replicate}.bam",
        bamindex="maplink/RNA-{condition}-{replicate}.bam.bai",
        annotation="annotation/ensembl.gtf"
    output:
        transcripts="spectre/transcripts/{condition}-{replicate}/transcripts.gtf",
        isoforms="spectre/transcripts/{condition}-{replicate}/isoforms.fpkm_tracking",
        genes="spectre/transcripts/{condition}-{replicate}/genes.fpkm_tracking",
    conda:
        "../envs/cufflinks.yaml"
    threads: 20
    shell:
        "mkdir -p spectre/transcripts; cufflinks {input.bam} -p {threads} -o ./spectre/transcripts/{wildcards.condition}-{wildcards.replicate}/ -g {input.annotation} --library-type fr-firststrand"

rule predictSpectre:
   input:
       bam="bam/RIBO-{condition}-{replicate}.bam",
       bamindex="maplink/RIBO-{condition}-{replicate}.bam.bai",
       gtf="annotation/ensembl.gtf",
       fpkm="spectre/transcripts/{condition}-{replicate}/isoforms.fpkm_tracking",
   output:
       results="spectre/{condition}-{replicate}/result.txt",
       log="spectre/{condition}-{replicate}/spectre_log.txt"
   threads: 5
   singularity:
        "docker://gelhausr/spectre:latest"
   shell:
       """
       mkdir -p spectre;
       SPECtre.py --input {input.bam} --output {output.results} --log {output.log} --fpkm {input.fpkm} --gtf {input.gtf} --nt {threads}
       """
