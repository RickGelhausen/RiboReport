def read_cutoff(filename):
    try:
        with open(filename, "r") as f:
            return f.readline().strip()
    except FileNotFoundError:
        return "failed"

rule createRibotricerIndex:
    input:
        genome=rules.retrieveGenome.output,
        annotation=rules.createEnsemblAnnotation.output
    output:
        "ribotricer/ribotricer_candidate_orfs.tsv"
    conda:
        "../envs/ribotricer.yaml"
    shell:
        """
        mkdir -p ribotricer;
        ribotricer prepare-orfs --gtf {input.annotation} --fasta {input.genome} --prefix ribotricer/ribotricer
        """

rule createRibotricerCutoff:
    input:
        bam_ribo="maplink/RIBO-{condition}-{replicate}.bam",
        ribo_index="maplink/RIBO-{condition}-{replicate}.bam.bai",
        bam_rna="maplink/RNA-{condition}-{replicate}.bam",
        rna_index="maplink/RNA-{condition}-{replicate}.bam.bai",
        index="ribotricer/ribotricer_candidate_orfs.tsv"
    output:
        "ribotricer/{condition}-{replicate}/recommended_cutoff.txt"
    conda:
        "../envs/ribotricer.yaml"
    shell:
        """
        mkdir -p ribotricer/{wildcards.condition}-{wildcards.replicate};
        ribotricer learn-cutoff --ribo_bams {input.bam_ribo} --rna_bams {input.bam_rna}  --prefix ribotricer/ribotricer --ribotricer_index {input.index} | tail -n 1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'> {output}
        """

rule predictRibotricer:
    input:
        bam_ribo="maplink/RIBO-{condition}-{replicate}.bam",
        index="ribotricer/ribotricer_candidate_orfs.tsv",
        cutoff="ribotricer/{condition}-{replicate}/recommended_cutoff.txt"
    output:
        "ribotricer/{condition}-{replicate}/ribotricer_translating_ORFs.tsv"
    conda:
        "../envs/ribotricer.yaml"
    params:
        cutoff=lambda wildcards, input: read_cutoff(input[2])
    shell:
        """
        mkdir -p ribotricer/{wildcards.condition}-{wildcards.replicate}
        ribotricer detect-orfs --bam {input.bam_ribo} --ribotricer_index {input.index} --prefix ribotricer/{wildcards.condition}-{wildcards.replicate}/ribotricer --phase_score_cutoff {params.cutoff}
        """

