rule fastqcunique:
    input:
        sam="sam/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/4unique/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/4unique/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/map/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "mkdir -p qc/4unique; fastqc -o qc/4unique -t {threads} -f sam_mapped {input.sam}; mv qc/4unique/{params.prefix}_fastqc.html {output.html}; mv qc/4unique/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcmulti:
    input:
        sam="sammulti/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/3mapped/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/3mapped/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/sammulti/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "mkdir -p qc/3mapped; fastqc -o qc/3mapped -t {threads} -f sam_mapped {input.sam}; mv qc/3mapped/{params.prefix}_fastqc.html {output.html}; mv qc/3mapped/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcraw:
    input:
        fastq="trimlink/{method}-{condition}-{replicate}.fastq.gz",
    output:
        html="qc/1raw/{method}-{condition}-{replicate}-raw_fastqc.html",
        zip="qc/1raw/{method}-{condition}-{replicate}-raw_fastqc.zip"
        #report("qc/1raw/{method}-{condition}-{replicate}-raw.html", caption="../report/fastqcraw.rst", category="Input quality control")
    conda:
        "../envs/fastqc.yaml"
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.fastq))[0])[0])
    threads: 8
    shell:
        "mkdir -p qc/1raw; fastqc -o qc/1raw -t {threads} {input.fastq}; mv qc/1raw/{params.prefix}_fastqc.html {output.html}; mv qc/1raw/{params.prefix}_fastqc.zip {output.zip}"

rule fastqctrimmed:
    input:
        reads="trimmed/{method}-{condition}-{replicate}.fastq"
    output:
        html="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html",
        zip="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.zip"
        #report("qc/2trimmed/{method}-{condition}-{replicate}-trimmed.html", caption="../report/fastqctrimmed.rst", category="Trimming")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "mkdir -p qc/2trimmed; fastqc -o qc/2trimmed -t {threads} {input}; mv qc/2trimmed/{params.prefix}_fastqc.html {output.html}; mv qc/2trimmed/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcrrnafilter:
    input:
        reads="bam/{method}-{condition}-{replicate}.bam"
    output:
        html="qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html",
        zip="qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA_fastqc.zip"
        #report("qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA.html", caption="../report/fastqcnorRNA.rst", category="Removing hits mapping to rRNA")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "mkdir -p qc/5removedrRNA; fastqc -o qc/5removedrRNA -t {threads} {input}; mv qc/5removedrRNA/{params.prefix}_fastqc.html {output.html}; mv qc/5removedrRNA/{params.prefix}_fastqc.zip {output.zip}"



rule featurescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/all/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        """
        mkdir -p qc/all;
        column3=$(cut -f3 auxiliary/unambigous_annotation.gff | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "gene" ]];
        then
            featureCounts -T {threads} -t gene -g ID -a {input.annotation} -o {output.txt} {input.bam};
        else
            touch {output.txt};
        fi
        """

rule trnafeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/trnainall/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        """
        mkdir -p qc/trnainall;
        column3=$(cut -f3 auxiliary/unambigous_annotation.gff | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "tRNA" ]];
        then
            featureCounts -T {threads} -t tRNA -g ID -a {input.annotation} -o {output.txt} {input.bam};
        else
            touch {output.txt};
        fi
        """

rule norrnafeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainall/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        """
        mkdir -p qc/rrnainall;
        column3=$(cut -f3 auxiliary/unambigous_annotation.gff | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "rRNA" ]];
        then
            featureCounts -T {threads} -t rRNA -g ID -a {input.annotation} -o {output.txt} {input.bam};
        else
            touch {output.txt};
        fi
        """

rule rrnatotalfeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bammulti/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainallaligned/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        """
        mkdir -p qc/rrnainallaligned;
        column3=$(cut -f3 auxiliary/unambigous_annotation.gff | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "rRNA" ]];
        then
            featureCounts -T {threads} -t rRNA -g ID -a {input.annotation} -o {output.txt} {input.bam};
        else
            touch {output.txt};
        fi
        """
        
rule rrnauniquefeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="rRNAbam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainuniquelyaligned/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        """
        mkdir -p qc/rrnainuniquelyaligned;
        column3=$(cut -f3 auxiliary/unambigous_annotation.gff | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "rRNA" ]];
        then
            featureCounts -T {threads} -t rRNA -g ID -a {input.annotation} -o {output.txt} {input.bam};
        else
            touch {output.txt};
        fi
        """

rule coveragedepth:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "coverage/{method}-{condition}-{replicate}.bed"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p coverage; bedtools genomecov -ibam {input} -bg > {output}"

rule multiqc:
    input:
        expand("qc/1raw/{method}-{condition}-{replicate}-raw_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/2trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/4unique/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/3mapped/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/all/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/trnainall/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainallaligned/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainuniquelyaligned/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainall/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("trimmed/{method}-{condition}-{replicate}.fastq", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
    output:
        report("qc/multi/multiqc_report.html", caption="../report/multiqc.rst", category="Quality control")
    params:
        dir="qc/multi"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "export LC_ALL=en_US.utf8; export LANG=en_US.utf8; multiqc -f -d --exclude picard --exclude gatk -z -o {params.dir} qc/3mapped qc/1raw qc/2trimmed qc/5removedrRNA qc/4unique qc/all qc/trnainall qc/rrnainallaligned qc/rrnainuniquelyaligned qc/rrnainall trimmed  2> {log}"
