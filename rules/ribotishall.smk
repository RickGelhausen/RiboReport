rule ribotishQualityRIBO:
    input:
        fp="maplink/RIBO/{condition}-{replicate}.bam",
        genome=rules.retrieveGenome.output,
        annotation=rules.ribotishAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        bamindex="maplink/RIBO/{condition}-{replicate}.bam.bai"
    output:
        reportpdf="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.pdf",
        reporttxt=report("ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.txt", caption="../report/ribotishquality.rst", category="Ribotish"),
        offsetdone="maplink/RIBO/{condition, [a-zA-Z]+}-{replicate,\d+}.qualdone"
    params:
        offsetparameters="maplink/RIBO/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.para.py"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotishquality.log"
    shell:
        "conda activate /beegfs/work/fr_fe1017/miniconda3/envs/ribotish2; mkdir -p ribotish; ribotish quality -v -p {threads} -b {input.fp} -g {input.annotation} -o {params.reporttxt} -f {params.reportpdf} 2> {log} || true; if grep -q \"offdict = {{'m0': {{}}}}\" {params.offsetparameters}; then mv {params.offsetparameters} {params.offsetparameters}.unused; fi; touch {output.offsetdone}"

rule ribotishQualityTIS:
    input:
        tis="maplink/TIS/{condition}-{replicate}.bam",
        genome=rules.retrieveGenome.output,
        annotation=rules.ribotishAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        tisindex="maplink/TIS/{condition}-{replicate}.bam.bai",
    output:
        offsetdone="maplink/TIS/{condition, [a-zA-Z]+}-{replicate,\d+}.qualdone"
    params:
        offsetparameters="maplink/TIS/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.para.py",
        reporttxt="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual_tis.txt",
        reportpdf="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual_tis.pdf"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotishqualitytis.log"
    shell:
        "conda activate /beegfs/work/fr_fe1017/miniconda3/envs/ribotish2; mkdir -p ribotish; ribotish quality -v -p {threads} -b {input.tis} -b {input.tis} -g {input.annotation} -o {params.reporttxt} -f {params.reportpdf} 2> {log} || true; if grep -q \"offdict = {{'m0': {{}}}}\" {params.offsetparameters}; then mv {params.offsetparameters} {params.offsetparameters}.unused; fi; touch {output.offsetdone}"

rule ribotish:
    input:
        fp= lambda wildcards: expand("maplink/RIBO/{{condition}}-{replicate}.bam", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        tis= lambda wildcards: expand("maplink/TIS/{{condition}}-{replicate}.bam", zip, replicate=samples.loc[(samples["method"] == "TIS") & (samples["condition"] == wildcards.condition), "replicate"]),
        genome=rules.retrieveGenome.output,
        annotation=rules.ribotishAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        bamindex= lambda wildcards: expand("maplink/RIBO/{{condition}}-{replicate}.bam.bai", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        tisindex= lambda wildcards: expand("maplink/TIS/{{condition}}-{replicate}.bam.bai", zip, replicate=samples.loc[(samples["method"] == "TIS") & (samples["condition"] == wildcards.condition), "replicate"]),
        #offsetparameters= lambda wildcards: expand("maplink/RIBO/{{condition}}-{replicate}.qualdone", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        report="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv_all.txt",
        #report=report("ribotish/{condition, [a-zA-Z]+}-newORFs.tsv_all.txt", caption="../report/ribotish.rst", category="Ribotish"),
        filtered="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv"
    params:
        fplist= lambda wildcards, input: ','.join(list(set(input.fp))),
        tislist= lambda wildcards, input: ','.join(list(set(input.tis))),
        codons= lambda wildcards: ("" if not CODONS else (" --alt --altcodons " + CODONS)),
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}_ribotish.log"
    shell:
        "conda activate /beegfs/work/fr_fe1017/miniconda3/envs/ribotish2; mkdir -p ribotish; ribotish predict --harr -v {params.codons} -p {threads} -t {params.tislist} -b {params.fplist} -g {input.annotation} -f {input.genome} -o {output.filtered} 2> {log}"
