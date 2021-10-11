
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
