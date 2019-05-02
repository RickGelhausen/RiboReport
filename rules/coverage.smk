rule coverage:
    input:
        bam="maplink/RIBO-{condition}-{replicate}.bam",
        bai="maplink/RIBO-{condition}-{replicate}.bam.bai"
    output:
        covfwd="coverage/{condition}-{replicate}_cov_fwd.bedgraph",
        covrev="coverage/{condition}-{replicate}_cov_rev.bedgraph",
        asitefwd="coverage/{condition}-{replicate}_asite_fwd.bedgraph",
        asiterev="coverage/{condition}-{replicate}_asite_rev.bedgraph"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p coverage; ribo_benchmark/scripts/coverage.py --alignment_file {input.bam} --output_file_prefix coverage/{wildcards.condition}-{wildcards.replicate}"
