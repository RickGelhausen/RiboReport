rule fastqcunique:
    input:
        bam="maplink/{method}-{condition}-{replicate}.bam"
    output:
        html="qc/unique/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/unique/{method}-{condition}-{replicate}-map_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.bam))[0])
    shell:
        "mkdir -p qc/unique; fastqc -o qc/unique -t {threads} -f bam_mapped {input.bam}; mv qc/unique/{params.prefix}_fastqc.html {output.html}; mv qc/unique/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcmapped:
    input:
        sam="sam/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/map/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/map/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/map/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "mkdir -p qc/map; fastqc -o qc/map -t {threads} -f sam_mapped {input.sam}; mv qc/map/{params.prefix}_fastqc.html {output.html}; mv qc/map/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcraw:
    input:
        reads=getfastq,
    output:
        html="qc/raw/{method}-{condition}-{replicate}-raw_fastqc.html",
        zip="qc/raw/{method}-{condition}-{replicate}-raw_fastqc.zip"
        #report("qc/raw/{method}-{condition}-{replicate}-raw.html", caption="../report/fastqcraw.rst", category="Input quality control")
    conda:
        "../envs/fastqc.yaml"
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.reads[0]))[0])[0])
    threads: 8
    shell:
        "mkdir -p qc/raw; fastqc -o qc/raw -t {threads} {input}; mv qc/raw/{params.prefix}_fastqc.html {output.html}; mv qc/raw/{params.prefix}_fastqc.zip {output.zip}"

rule fastqctrimmed:
    input:
        reads="trimmed/{method}-{condition}-{replicate}.fastq"
    output:
        html="qc/trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html",
        zip="qc/trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.zip"
        #report("qc/trimmed/{method}-{condition}-{replicate}-trimmed.html", caption="../report/fastqctrimmed.rst", category="Trimming")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "mkdir -p qc/trimmed; fastqc -o qc/trimmed -t {threads} {input}; mv qc/trimmed/{params.prefix}_fastqc.html {output.html}; mv qc/trimmed/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcrrnafilter:
    input:
        reads="norRNA/{method}-{condition}-{replicate}.fastq"
    output:
        html="qc/norRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html",
        zip="qc/norRNA/{method}-{condition}-{replicate}-norRNA_fastqc.zip"
        #report("qc/norRNA/{method}-{condition}-{replicate}-norRNA.html", caption="../report/fastqcnorRNA.rst", category="Removing hits mapping to rRNA")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "mkdir -p qc/norRNA; fastqc -o qc/norRNA -t {threads} {input}; mv qc/norRNA/{params.prefix}_fastqc.html {output.html}; mv qc/norRNA/{params.prefix}_fastqc.zip {output.zip}"

rule multiqc:
    input:
        expand("qc/raw/{method}-{condition}-{replicate}-raw_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/norRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/map/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/unique/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
    output:
        report("qc/multi/multiqc_report.html", caption="../report/multiqc.rst", category="Quality control")
    params:
        dir="qc/multi"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "export LC_ALL=en_US.utf8; export LANG=en_US.utf8; multiqc -f -d --exclude picard --exclude gatk -z -o {params.dir} qc/raw qc/trimmed qc/norRNA qc/map qc/unique trimmed  2> {log}"
