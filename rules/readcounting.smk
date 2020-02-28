rule generateAnnotationReadCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="auxiliary/unambigous_annotation.gtf"
    output:
        "auxiliary/annotation_reads.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p auxiliary
        UNIQUE="$(cut -f3 {input.annotation} | sort | uniq)"
        IDENTIFIER="ID"
        LINE="$(sed '3q;d' {input.annotation})"
        if [[ $LINE == *"gene_id="* ]]; then IDENTIFIER="gene_id"; fi;
        for f in ${{UNIQUE}}
        do
            featureCounts -F GTF -s 1 -g $IDENTIFIER -O -t $f -M --fraction -a {input.annotation} {input.bam} -T {threads} -o auxiliary/annotation_total_reads.raw.tmp
            cat auxiliary/annotation_total_reads.raw.tmp | sed 1,2d | awk -v var=$f -FS'\\t' '{{print $0"\\t"var}}' >> {output}
            rm auxiliary/annotation_total_reads.raw.tmp
        done
        """

rule mapReads:
    input:
        reads="auxiliary/annotation_reads.raw",
        annotation="auxiliary/enriched_annotation.gtf"
    output:
        "auxiliary/total_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary; RiboReport/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule totalMappedReads:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        mapped="auxiliary/total_mapped_reads.txt",
        length="auxiliary/total_read_lengths.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; RiboReport/scripts/total_mapped_reads.py -b {input.bam} -m {output.mapped} -l {output.length}"

rule createExcelSummary:
    input:
        total="auxiliary/total_mapped_reads.txt",
        reads="auxiliary/total_annotation.gtf",
        genome="genomes/genome.fa"
    output:
        xlsx="auxiliary/final_annotation.xlsx",
        annotation="auxiliary/final_annotation.gff",
        complete="auxiliary/final_annotation_complete.gff"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; RiboReport/scripts/generate_output_annotation.py -t {input.total} -r {input.reads} -g {input.genome} -o {output.xlsx} --simple {output.annotation} --complete {output.complete}"
