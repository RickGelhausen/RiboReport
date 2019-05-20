rule reparationGFF:
    input:
        "reparation/{condition}-{replicate}/Predicted_ORFs.txt"
    output:
        "reparation/{condition, [a-zA-Z]+}-{replicate,\d+}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; ribo_benchmark/scripts/reparationGFF.py -c {wildcards.condition}  -i {input} -o {output}"

rule concatReparation:
    input:
        lambda wildcards: expand("reparation/{{condition}}-{replicate}.reparation.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; ribo_benchmark/scripts/concatGFF.py {input} -o {output}"

rule concatIrsom:
    input:
        lambda wildcards: expand("irsom/{{condition}}-{replicate}.irsom.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.irsom.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; ribo_benchmark/scripts/concatGFF.py {input} -o {output}"

rule deepriboGFF:
    input:
        "deepribo/{condition}-{replicate}/predictions.csv"
    output:
        "deepribo/{condition, [a-zA-Z]+}-{replicate,\d+}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; ribo_benchmark/scripts/deepriboGFF.py -c {wildcards.condition}  -i {input} -o {output}"

rule concatDeepRibo:
    input:
        lambda wildcards: expand("deepribo/{{condition}}-{replicate}.deepribo.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; ribo_benchmark/scripts/concatGFF.py {input} -o {output}"

rule ribotishGFF:
    input:
        "ribotish/{condition}-newORFs.tsv_all.txt"
    output:
        "tracks/{condition, [a-zA-Z]+}.ribotish.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; ribo_benchmark/scripts/ribotish.py {input} --condition {wildcards.condition} --output_gff3_filepath {output}"

rule mergeConditions:
    input:
        ribotish="tracks/{condition}.ribotish.gff",
        reparation="tracks/{condition}.reparation.gff",
        deepribo="tracks/{condition}.deepribo.gff",
        irsom="tracks/{condition}.irsom.gff"
    output:
        report("tracks/{condition, [a-zA-Z]+}.merged.gff", caption="../report/novelmerged.rst", category="Novel ORFs")
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; cat {input.ribotish} > {output}.unsorted; cat {input.reparation} >> {output}.unsorted; cat {input.deepribo} >> {output}.unsorted; cat {input.irsom} >> {output}.unsorted; bedtools sort -i {output}.unsorted > {output};"

rule mergeAll:
    input:
        mergedGff=expand("tracks/{condition}.merged.gff", zip, condition=set(samples["condition"]))
    output:
        report("tracks/all.gff", caption="../report/novelall.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; ribo_benchmark/scripts/concatGFF.py {input.mergedGff} -o {output}"

rule filterAll:
    input:
        "tracks/all.gff"
    output:
        report("tracks/combined.gtf", caption="../report/combined.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; ribo_benchmark/scripts/allToUniqueGTF.py -i {input} -o {output}"


rule intersectOperonRef:
    input:
        reference="gtfdata/all_filtered_massspec.gtf",
        operons="gtfdata/Ecoli_list_of_operons_V1.gtf"
    output:
        "intersect/ref_operon_intersect.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -a {input.reference} -b {input.operons} -u -s > {output}"


rule intersectOperonTools:
    input:
        tools="tracks/combined.gtf",
        operons="gtfdata/Ecoli_list_of_operons_V1.gtf"
    output:
        "intersect/tools_operon_intersect.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -a {input.tools} -b {input.operons} -u -s > {output}"

rule intersectOperonRefComplement:
    input:
        reference="gtfdata/all_filtered_massspec.gtf",
        operons="gtfdata/Ecoli_list_of_operons_V1.gtf"
    output:
        "intersect/ref_operon_complement.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -v -a {input.reference} -b {input.operons} -s > {output}"

rule intersectOperonToolsComplement:
    input:
        tools="tracks/combined.gtf",
        operons="gtfdata/Ecoli_list_of_operons_V1.gtf"
    output:
        "intersect/tools_operon_complement.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -v -a {input.tools} -b {input.operons} -s > {output}"

rule intersectSmallRef:
    input:
        reference="gtfdata/all_filtered_massspec.gtf",
        small="gtfdata/smallORFs.gtf"
    output:
        "intersect/ref_smallORFs_intersect.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -a {input.reference} -b {input.small} -u -s > {output}"

rule intersectSmallTools:
    input:
        tools="tracks/combined.gtf",
        small="gtfdata/smallORFs.gtf"
    output:
        "intersect/tools_smallORFs_intersect.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -a {input.tools} -b {input.small} -u -s > {output}"

rule intersectncRNARef:
    input:
        reference="gtfdata/all_filtered_massspec.gtf",
        small="gtfdata/ncRNAs.gtf"
    output:
        "intersect/ref_ncRNAs_intersect.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -a {input.reference} -b {input.small} -u -s > {output}"

rule intersectncRNATools:
    input:
        tools="tracks/combined.gtf",
        small="gtfdata/ncRNAs.gtf"
    output:
        "intersect/tools_ncRNAs_intersect.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -a {input.tools} -b {input.small} -u -s > {output}"

rule intersectPseudoGenesRef:
    input:
        reference="gtfdata/all_filtered_massspec.gtf",
        pseudo="gtfdata/pseudogenes.gtf"
    output:
        "intersect/ref_pseudogenes_intersect.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -a {input.reference} -b {input.pseudo} -u -s > {output}"

rule intersectPseudoGenesTools:
    input:
        tools="tracks/combined.gtf",
        pseudo="gtfdata/pseudogenes.gtf"
    output:
        "intersect/tools_pseudogenes_intersect.gtf"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p intersect; bedtools intersect -a {input.tools} -b {input.pseudo} -u -s > {output}"

