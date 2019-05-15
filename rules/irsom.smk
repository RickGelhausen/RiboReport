rule trainIrsom:
    input:
        coding="ecoli/coding.fa",
        noncoding="ecoli/noncoding.fa"
    output:
        "irsom/model/test.txt"
    threads: 1
    shell:
        """
        source activate irsom;
        mkdir -p irsom;
        python tools/IRSOM/scripts/train.py --featurer=tools/IRSOM/bin/Featurer -c {input.coding} -n {input.noncoding} --output=irsom/model
        """
