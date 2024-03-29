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

rule smorferFourierTransformWithCalibratedBam:
    input:
        bam="calibrated_bam/RIBO-{condition}-{replicate}_calibrated.bam",
        bed="smorfer/{condition}-{replicate}/RPF_high.bed"
    output:
        "smorfer/{condition}-{replicate}/RPF_3nt_translated.txt"
    singularity:
        "docker://gelhausr/smorfer:latest"
    threads: 1
    shell:
        """
        mkdir -p smorfer/{wildcards.condition}-{wildcards.replicate};
        smorfer.sh -s FT_RPF.sh high_expressed_ORFs.bed calibrated_RiboSeq.bam -b {input.bed} -a {input.bam} -o smorfer/{wildcards.condition}-{wildcards.replicate}/
        """

rule smorferFilterStartWithCalibratedBam:
    input:
        bam="calibrated_bam/RIBO-{condition}-{replicate}_calibrated.bam",
        bed="smorfer/{condition}-{replicate}/RPF_translated.txt"
    output:
        "smorfer/{condition}-{replicate}/best_start_results.txt"
    singularity:
        "docker://gelhausr/smorfer:latest"
    threads: 1
    shell:
        """
        mkdir -p smorfer/{wildcards.condition}-{wildcards.replicate};
        smorfer.sh -s find_best_start.sh -b {input.bed} -a {input.bam} -o smorfer/{wildcards.condition}-{wildcards.replicate}/
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

rule smorferGFF:
    input:
        #"smorfer/{condition}-{replicate}/best_start_results.txt"
        "smorfer/{condition}-{replicate}/RPF_translated.txt"
    output:
        "smorfer/{condition, [a-zA-Z]+}-{replicate,\d+}.smorfer.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/smorferGFF.py -c {wildcards.condition}  -i {input} -o {output}"

rule concatSmorfer:
    input:
        lambda wildcards: expand("smorfer/{{condition}}-{replicate}.smorfer.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.smorfer.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/concatGFF.py {input} -o {output}"