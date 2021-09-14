def read_first_chrom(filename):
    try:
        line = ""
        with open(filename, "r") as f:
            line = f.readline().split()[0][1:]
        return line
    except FileNotFoundError:
        return "failed"

rule createSmorferPutativeORFs:
    input:
        genome=rules.retrieveGenome.output,
    output:
        "smorfer/smorfer_pORFs.bed"
    singularity:
        "docker://gelhausr/smorfer:latest"
    threads: 1
    params:
        genome_name=lambda wildcards, input: read_first_chrom(input[0])
    shell:
        """
        mkdir -p smorfer;
        smorfer.sh -s putative_orfs.sh -g {input.genome} -c {params.genome_name} -n smorfer_pORFs -i 9-150 -o smorfer/
        """

rule createSmorferFourierTransform:
    input:
        genome=rules.retrieveGenome.output,
        bed="smorfer/smorfer_pORFs.bed"
    output:
        "smorfer/FT_passed.bed"
    singularity:
        "docker://gelhausr/smorfer:latest"
    threads: 1
    shell:
        """
        mkdir -p smorfer;
        smorfer.sh -s FT_GCcontent.sh -g {input.genome} -o smorfer/ -b {input.bed}
        """

rule createSmorferRPFcount:
    input:
        bam="maplink/RIBO-{condition}-{replicate}.bam",
        bed="smorfer/FT_passed.bed"
    output:
        "smorfer/{condition}-{replicate}/RPF_translated.txt"
    singularity:
        "docker://gelhausr/smorfer:latest"
    threads: 1
    shell:
        """
        mkdir -p smorfer/{wildcards.condition}-{wildcards.replicate};
        smorfer.sh -s count_RPF.sh -b {input.bed} -a {input.bam} -o smorfer/{wildcards.condition}-{wildcards.replicate}/
        """


rule createSmorferRPFcountNOFT:
    input:
        bam="maplink/RIBO-{condition}-{replicate}.bam",
        bed="smorfer/smorfer_pORFs.bed"
    output:
        "smorfer/{condition}-{replicate}_noFT/RPF_translated.txt"
    singularity:
        "docker://gelhausr/smorfer:latest"
    threads: 1
    shell:
        """
        mkdir -p smorfer/{wildcards.condition}-{wildcards.replicate}_noFT;
        smorfer.sh -s count_RPF.sh -b {input.bed} -a {input.bam} -o smorfer/{wildcards.condition}-{wildcards.replicate}_noFT/
        """
