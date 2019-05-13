rule createCodingBed:
    input:
        "annotation/annotation.gtf"
    output:
        coding="cpat/annotation_coding.bed",
        noncoding="cpat/annotation_noncoding.bed"
        all="cpat/annotation_all.bed"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p cpat; ribo_benchmark/scripts/retrieveCodingInformation.py -i {input} -o cpat/annotation"

rule bedToFasta:
    input:
        genome="genomes/genome.fa",
        coding="cpat/annotation_coding.bed",
        noncoding="cpat/annotation_noncoding.bed"
    output:
        coding="cpat/annotation_coding.fa",
        noncoding"cpat/annotation_noncoding.fa"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        """
        mkdir -p cpat;
        bedtools getfasta -fi {input.genome} -bed {input.coding} -fo {output.coding}
        bedtools getfasta -fi {input.genome} -bed {input.noncoding} -fo {output.noncoding}
        """

rule makeHexamerTab:
    input:
        coding="cpat/annotation_coding.fa",
        noncoding="cpat/annotation_noncoding.fa"
    output:
        "cpat/hexamer.tsv"
    threads: 1
    shell:
        """
        conda activate cpat;
        mkdir -p cpat;
        make_hexamer_tab.py -c {input.coding} -n {input.noncoding} > {output};
        conda deactivate cpat;
        """

rule makeLogitModel:
    input:
        hexamer="cpat/hexamer.tsv",
        coding="cpat/annotation_coding.fa",
        noncoding="cpat/annotation_noncoding.fa"
    output:
        xls="cpat/model.feature.xls",
        rdata="cpat/model.logit.RData",
        model="cpat/model.make_logitModel.r"
    threads: 1
    shell:
        """
        conda activate cpat;
        make_logitModel.py -x {input.hexamer} -c {input.coding} -n {input.noncoding} -o cpat/model;
        conda deactivate cpat;
        """

rule cpat:
    input:
        genome="genomes/genome.fa",
        annotation="cpat/annotation_all.bed",
        hexamer="cpat/hexamer.tsv",
        logit="cpat/model.logit.RData"
    output:
        "cpat/cpat_predictions.tsv"
    threads: 1
    shell:
        """
        conda activate cpat;
        cpat.py -r {input.genome} -g {input.annotation} -x {input.hexamer} -d {input.logit} -o {output}
        conda deactivate cpat;
        """
