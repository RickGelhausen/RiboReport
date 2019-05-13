rule generateMetageneRoi:
    input:
        annotation="annotation/annotation_processed.gtf"
    output:
        "offsets/metagene_rois.txt"
    conda:
        "../envs/auxiliary.yaml"
    threads: 1
    shell:
        "mkdir -p offsets; metagene generate -q offsets/metagene --landmark cds_start --annotation_files {input.annotation}"

rule psiteOffsets:
    input:
        rois=rules.generateMetageneRoi.output,
        bam="maplink/{method}-{condition}-{replicate}.bam",
        bamindex="maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        psite="offsets/{method}-{condition}-{replicate}_p_offsets.txt"
    conda:
        "../envs/auxiliary.yaml"
    threads: 1
    shell:
        "mkdir -p offsets; psite -q {input.rois} offsets/{wildcards.method}-{wildcards.condition}-{wildcards.replicate} --min_length 22 --max_length 40 --require_upstream --count_files {input.bam}; sed -i '/^#/ d' {output}"

rule gff2gtf:
    input:
        annotation={rules.retrieveAnnotation.output},
    output:
        gtf="auxiliary/featurecount/annotation.gtf",
    conda:
        "../envs/cufflinks.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary/featurecount; gffread {input.annotation} -T -o {output.gtf}"
