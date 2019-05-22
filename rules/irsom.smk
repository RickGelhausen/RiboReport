rule trainIrsom:
    input:
        coding="ecoli/coding.fasta",
        noncoding="ecoli/noncoding.fasta"
    output:
        "irsom/model/test.txt"
    threads: 1
    shell:
        """
        source activate irsom
        mkdir -p irsom;
        python tools/IRSOM/scripts/train.py --featurer=tools/IRSOM/bin/Featurer -c {input.coding} -n {input.noncoding} --output=irsom/model
        """

rule bedtoolsGetfasta:
    input:
        genome="genomes/genome.fa",
        transcripts="transcripts/{condition}-{replicate}/transcripts.bed"
    output:
        "transcripts/{condition}-{replicate}/transcripts.fa"
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    shell:
        "mkdir -p transcripts; bedtools getfasta -fi {input.genome} -name -bed {input.transcripts} -fo {output}"

rule predictIrsom:
    input:
        transcripts="transcripts/{condition}-{replicate}/transcripts.fa",
    output:
        "irsom/{condition}-{replicate}/test.txt"
    threads: 1
    shell:
        """
        source activate irsom
        mkdir -p irsom;
        python tools/IRSOM/scripts/predict.py --featurer=tools/IRSOM/bin/Featurer --file={input.transcripts} --model=tools/IRSOM/model/Escherichia_coli/ --output=irsom/{wildcards.condition}-{wildcards.replicate}/
        """
