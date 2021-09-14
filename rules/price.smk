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
