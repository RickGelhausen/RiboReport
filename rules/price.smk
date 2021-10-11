rule createPriceIndex:
    input:
        genome=rules.retrieveGenome.output,
        annotation=rules.createEnsemblAnnotation.output
    output:
        "price/annotation_ensembl.oml"
    singularity:
        "docker://gelhausr/price:latest"
    threads: 1
    shell:
        """
        mkdir -p price;
        gedi -e IndexGenome -s {input.genome} -a {input.annotation} -nobowtie -nokallisto -nostar -f price/ -o {output} -n price_genome
        """

rule predictPrice:
    input:
        genomic="price/annotation_ensembl.oml",
        bam="maplink/RIBO-{condition}-{replicate}.bam",
        bamindex="maplink/RIBO-{condition}-{replicate}.bam.bai"
    output:
        cit="price/{condition}-{replicate}/results.orfs.cit",
        tsv="price/{condition}-{replicate}/results.orfs.tsv"
    singularity:
        "docker://gelhausr/price:latest"
    threads: 1
    shell:
        """
        mkdir -p price;
        gedi -e Price -reads {input.bam} -genomic {input.genomic} -prefix price/{wildcards.condition}-{wildcards.replicate}/results
        """

rule decodeResultsPrice:
    input:
        cit="price/{condition}-{replicate}/results.orfs.cit"
    output:
        bed="price/{condition}-{replicate}/results.orfs.filtered.bed"
    singularity:
        "docker://gelhausr/price:latest"
    threads: 1
    shell:
        """
        mkdir -p price;
        gedi -e ViewCIT -m bed -name 'd.getType()' {input.cit} > {output.bed}
        """

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

