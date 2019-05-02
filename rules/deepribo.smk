rule parseDeepRibo:
    input:
        covS= "coverage/{method}-{condition}-{replicate}_cov_fwd.bedgraph",
        covAS= "coverage/{method}-{condition}-{replicate}_cov_rev.bedgraph",
        asiteS= "coverage/{method}-{condition}-{replicate}_asite_fwd.bedgraph",
        asiteAS= "coverage/{method}-{condition}-{replicate}_asite_rev.bedgraph",
        genome= rules.retrieveGenome.output,
        annotation= rules.retrieveAnnotation.output
    output:
        "deepribo/{method}-{condition}-{replicate}/data_list.csv"
    conda:
        "../envs/deepribo.yaml"
    threads: 1
    shell:
        """
        mkdir -p deepribo/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}/0/;
        mkdir -p deepribo/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}/1/;
        python tools/DeepRibo/src/DataParser.py {input.covS} {input.covAS} {input.asiteS} {input.asiteAS} {input.genome} deepribo/{wildcards.method}-{wildcards.condition}-{wildcards.replicate} -g {input.annotation}
        """

#rule trainDeepRibo:
#    input:
#        "deepribo/{method}-{condition}-{replicate}/data_list.csv"
#    output:
#        "deepribo/{method}-{condition}-{replicate}/pred/bla.txt"
#    conda:
#        "../envs/deepribo.yaml"
#    threads: 20
#    shell:
#        "mkdir -p deepribo; python tools/DeepRibo/src/DeepRibo.py train deepribo/ -r 0.12 -c 0.13 --train_data {wildcards.method}-{wildcards.condition}-{wildcards.replicate}/  --num_workers {threads}"

predictDeepRibo:
    input:
        model="tools/DeepRibo/article/model_coli.pt",
        data="deepribo/{method}-{condition}-{replicate}/data_list.csv"
    output:
        "deepribo/{method}-{condition}-{replicate}_predictions.csv"
    conda:
        "../envs/deepribo.yaml"
    threads: 20
    shell:
        "mkdir -p deepribo; python tools/DeepRibo/src/DeepRibo.py predict deepribo/ --pred_data {wildcards.method}-{wildcards.condition}-{wildcards.replicate}/ -r 0.23 -c 0.10 --num_workers {threads} --dest {output}"
