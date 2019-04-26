rule parseDeepRibo:
    input:
        covS= "coverage/{method}-{condition}-{replicate}_cov_fwd.bedgraph",
        covAS= "coverage/{method}-{condition}-{replicate}_cov_rev.bedgraph",
        asiteS= "coverage/{method}-{condition}-{replicate}_asite_fwd.bedgraph",
        asiteAS= "coverage/{method}-{condition}-{replicate}_asite_rev.bedgraph",
        genome= rules.retrieveGenome.output,
        annotation= rules.retrieveAnnotation.output
    output:
        "deepribo/{method}-{condition}-{replicate}/0/noclue"
    conda:
        "../envs/deepribo.yaml"
    threads: 1
    log:
        "logs/deepribo_dataparser.log"
    shell:
        """
        mkdir -p deepribo/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}/0/;
        mkdir -p deepribo/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}/1/;
        python tools/DeepRibo/src/DataParser.py {input.covS} {input.covAS} {input.asiteS} {input.asiteAS} {input.genome} deepribo/{wildcards.method}-{wildcards.condition}-{wildcards.replicate} -g {input.annotation}
        """
