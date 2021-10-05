
rule reparationGFF:
    input:
        "reparation/{condition}-{replicate}/Predicted_ORFs.txt"
    output:
        "reparation/{condition, [a-zA-Z]+}-{replicate,\d+}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/reparationGFF.py -c {wildcards.condition}  -r {wildcards.replicate} -i {input} -o {output}"

rule concatReparation:
    input:
        lambda wildcards: expand("reparation/{{condition}}-{replicate}.reparation.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/concatGFF.py {input} -o {output}"

rule deepriboGFF:
    input:
        "deepribo/{condition}-{replicate}/predictions.csv"
    output:
        "deepribo/{condition, [a-zA-Z]+}-{replicate,\d+}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/deepriboGFF.py -c {wildcards.condition}  -i {input} -o {output}"

rule concatDeepRibo:
    input:
        lambda wildcards: expand("deepribo/{{condition}}-{replicate}.deepribo.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/concatGFF.py {input} -o {output}"

rule irsomGFF:
    input:
        "irsom/{condition}-{replicate}/result.txt"
    output:
        "irsom/{condition, [a-zA-Z]+}-{replicate,\d+}.irsom.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/irsomGFF.py -c {wildcards.condition}  -i {input} -o {output}"

rule concatIrsom:
    input:
        lambda wildcards: expand("irsom/{{condition}}-{replicate}.irsom.gff", zip, replicate=samples.loc[(samples["method"] == "RNA") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.irsom.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/concatGFF.py {input} -o {output}"

rule spectreGFF:
    input:
        "spectre/{condition}-{replicate}/result.txt"
    output:
        "spectre/{condition, [a-zA-Z]+}-{replicate,\d+}.spectre.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/spectreGFF.py -c {wildcards.condition}  -i {input} -o {output}"

rule concatSpectre:
    input:
        lambda wildcards: expand("spectre/{{condition}}-{replicate}.spectre.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.spectre.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/concatGFF.py {input} -o {output}"

rule ribotricerGFF:
    input:
        "ribotricer/{condition}-{replicate}/ribotricer_translating_ORFs.tsv"
    output:
        "ribotricer/{condition, [a-zA-Z]+}-{replicate,\d+}.ribotricer.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/ribotricerGFF.py -c {wildcards.condition}  -i {input} -o {output}"

rule concatRibotricer:
    input:
        lambda wildcards: expand("ribotricer/{{condition}}-{replicate}.ribotricer.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.ribotricer.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/concatGFF.py {input} -o {output}"

rule smorferGFF:
    input:
        "smorfer/{condition}-{replicate}/best_start_results.txt"
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

rule priceGFF:
    input:
        filtered="price/{condition}-{replicate}/results.orfs.filtered.bed",
        unfiltered="price/{condition}-{replicate}/results.orfs.tsv"
    output:
        filtered="price/{condition, [a-zA-Z]+}-{replicate,\d+}.price.filtered.gff",
        unfiltered="price/{condition, [a-zA-Z]+}-{replicate,\d+}.price.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/priceGFF.py -c {wildcards.condition}  -i {input.filtered} -u {input.unfiltered} -o {output.unfiltered} -f {output.filtered}"

rule concatPrice:
    input:
        unfiltered=lambda wildcards: expand("price/{{condition}}-{replicate}.price.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        filtered=lambda wildcards: expand("price/{{condition}}-{replicate}.price.filtered.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        unfiltered="tracks/{condition, [a-zA-Z]+}.price.gff",
        filtered="tracks/{condition, [a-zA-Z]+}.price.filtered.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p tracks;
        RiboReport/scripts/concatGFF.py {input.unfiltered} -o {output.unfiltered}
        RiboReport/scripts/concatGFF.py {input.filtered} -o {output.filtered}
        """

rule ribotishGFF:
    input:
        "ribotish/{condition}-newORFs.tsv_all.txt"
    output:
        "tracks/{condition, [a-zA-Z]+}.ribotish.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/ribotishGFF.py -c {wildcards.condition} -i {input} -o {output}"

rule mergeConditions:
    input:
        ribotish="tracks/{condition}.ribotish.gff",
        reparation="tracks/{condition}.reparation.gff",
        deepribo="tracks/{condition}.deepribo.gff",
        irsom="tracks/{condition}.irsom.gff",
        spectre="tracks/{condition}.spectre.gff",
        smorfer="tracks/{condition}.smorfer.gff",
        ribotricer="tracks/{condition}.ribotricer.gff",
        price="tracks/{condition}.price.gff"
    output:
        "tracks/{condition, [a-zA-Z]+}.merged.gff"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        """
        mkdir -p tracks;
        cat {input.spectre} > {output}.unsorted;
        cat {input.ribotish} >> {output}.unsorted;
        cat {input.reparation} >> {output}.unsorted;
        cat {input.deepribo} >> {output}.unsorted;
        cat {input.irsom} >> {output}.unsorted;
        cat {input.smorfer} >> {output}.unsorted;
        cat {input.ribotricer} >> {output}.unsorted;
        cat {input.price} >> {output}.unsorted;
        bedtools sort -i {output}.unsorted > {output};
        """

rule mergeAll:
    input:
        mergedGff=expand("tracks/{condition}.merged.gff", zip, condition=set(samples["condition"]))
    output:
        "tracks/all.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/concatGFF.py {input.mergedGff} -o {output}"

rule filterAll:
    input:
        "tracks/all.gff"
    output:
        "tracks/predictions.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/allToUniqueGTF.py -i {input} -o {output}"
