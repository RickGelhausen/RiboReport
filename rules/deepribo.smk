rule parseDeepRibo:
    input:
        covS= "coverage/{condition}-{replicate}_cov_fwd.bedgraph",
        covAS= "coverage/{condition}-{replicate}_cov_rev.bedgraph",
        asiteS= "coverage/{condition}-{replicate}_asite_fwd.bedgraph",
        asiteAS= "coverage/{condition}-{replicate}_asite_rev.bedgraph",
        genome= rules.retrieveGenome.output,
        annotation= rules.retrieveAnnotation.output
    output:
        "deepribo/{condition}-{replicate}/data_list.csv"
    conda:
        "../envs/deepribo.yaml"
    threads: 1
    shell:
        """
        mkdir -p deepribo/{wildcards.condition}-{wildcards.replicate}/0/;
        mkdir -p deepribo/{wildcards.condition}-{wildcards.replicate}/1/;
        python tools/DeepRibo/src/DataParser.py {input.covS} {input.covAS} {input.asiteS} {input.asiteAS} {input.genome} deepribo/{wildcards.condition}-{wildcards.replicate} -g {input.annotation}
        """
        
rule parameterEstimation:
    input:
        "deepribo/{condition}-{replicate}/data_list.csv"
    output:
        "deepribo/{condition}-{replicate}/figure/dunno.txt"
    conda:
        "../envs/estimation.yaml"
    threads: 1
    shell:
        "mkdir -p deepribo; Rscript ribo_benchmark/scripts/parameter_estimation.R -f {input}"


rule predictDeepRibo:
    input:
        model= "tools/DeepRibo/models/DeepRibo_model_v1.pt",
        data= "deepribo/{condition}-{replicate}/data_list.csv"
    output:
        "deepribo/{condition}-{replicate}/predictions.csv"
    conda:
        "../envs/deepribo.yaml"
    threads: 10
    shell:
        """
        mkdir -p deepribo;
        python tools/DeepRibo/src/DeepRibo.py predict deepribo/ --pred_data {wildcards.condition}-{wildcards.replicate}/ -r 0.2440181 -c 0.127214 --model {input.model} --dest {output} --num_workers {threads}
        """
