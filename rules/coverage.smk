rule coverage:
    input:
        bam="maplink/{method}-{condition}-{replicate}.bam",
        bai="maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        covfwd="coverage/{method}-{condition}-{replicate}_cov_fwd.bedgraph",
        covrev="coverage/{method}-{condition}-{replicate}_cov_rev.bedgraph",
        asitefwd="coverage/{method}-{condition}-{replicate}_asite_fwd.bedgraph",
        asiterev="coverage/{method}-{condition}-{replicate}_asite_rev.bedgraph"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p coverage; ribo_benchmark/scripts/coverage.py --alignment_file {input.bam} --output_file_prefix coverage/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}"
