rule coverageDepthFwd:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "coverage/{method}-{condition}-{replicate}_cov_fwd.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p coverage; bedtools genomecov -ibam {input} -d -strand + > {output}"

rule coverageDepthRev:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "coverage/{method}-{condition}-{replicate}_cov_rev.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p coverage; bedtools genomecov -ibam {input} -d -strand - > {output}"

rule aSiteOccupancy:
    input:
        bam="maplink/{method}-{condition}-{replicate}.bam",
        bai="maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        fwd="coverage/{method}-{condition}-{replicate}_asite_fwd.bedgraph",
        rev="coverage/{method}-{condition}-{replicate}_asite_rev.bedgraph"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p coverage; ribo_benchmark/scripts/generateASiteOccupancy.py --alignment_file {input.bam} --output_file_prefix coverage/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}"
