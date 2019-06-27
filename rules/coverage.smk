rule coverage:
    input:
        bam="maplink/RIBO-{condition}-{replicate}.bam",
        bai="maplink/RIBO-{condition}-{replicate}.bam.bai"
    output:
        covfwd="coverage_deepribo/{condition}-{replicate}_cov_fwd.bedgraph",
        covrev="coverage_deepribo/{condition}-{replicate}_cov_rev.bedgraph",
        asitefwd="coverage_deepribo/{condition}-{replicate}_asite_fwd.bedgraph",
        asiterev="coverage_deepribo/{condition}-{replicate}_asite_rev.bedgraph"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p coverage_deepribo; ribo_benchmark/scripts/coverage_deepribo.py --alignment_file {input.bam} --output_file_prefix coverage_deepribo/{wildcards.condition}-{wildcards.replicate}"
