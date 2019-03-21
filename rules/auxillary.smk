rule processAnnotation:
    input:
        annotation=rules.retrieveAnnotation.output
    output:
        "annotation/processed-annotation.gtf"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p annotation; python3 scripts/processAnnotation.py -a {input.annotation} -o {output}"

rule generateMetageneRoi:
    input:
        annotation=rules.processAnnotation.output
    output:
        "plastid/metagene_rois.txt"
    conda:
        "../envs/auxillary.yaml"
    threads: 1
    shell:
        "mkdir -p offsets; metagene generate -q offsets/metagene --landmark cds_start --annotation_files {input.annotation}"

rule psiteOffsets:
    input:
        rois=rules.generateMetageneRoi.output,
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        psite="offsets/{method}-{condition}-{replicate}_p_offsets.txt"
    conda:
        "../envs/auxillary.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p offsets; psite {input.rois} offsets/{params.prefix} --min_length 22 --max_length 40 --require_upstream --count_files {input.bam}"
