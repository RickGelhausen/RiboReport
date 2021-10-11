
rule transcriptsBED:
    input:
        "transcripts/{condition}-{replicate}/transcripts.gtf"
    output:
        "transcripts/{condition}-{replicate}/transcripts.bed"
    threads: 1
    conda:
        "../envs/mergetools.yaml"
    shell:
        "mkdir -p transcripts; RiboReport/scripts/transcriptsToBed.py -i {input} -o {output}"

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

rule irsomGetModel:
    input:
        slsomw=HTTP.remote("https://forge.ibisc.univ-evry.fr/lplaton/IRSOM/raw/master/model/Escherichia_coli/SLSOM/W.txt", keep_local=True),
        slsomb=HTTP.remote("https://forge.ibisc.univ-evry.fr/lplaton/IRSOM/raw/master/model/Escherichia_coli/SLSOM/biases.txt", keep_local=True),
        slsomp=HTTP.remote("https://forge.ibisc.univ-evry.fr/lplaton/IRSOM/raw/master/model/Escherichia_coli/SLSOM/parameters.txt", keep_local=True),
        somu=HTTP.remote("https://forge.ibisc.univ-evry.fr/lplaton/IRSOM/raw/master/model/Escherichia_coli/SOM/units.txt", keep_local=True),
        somp=HTTP.remote("https://forge.ibisc.univ-evry.fr/lplaton/IRSOM/raw/master/model/Escherichia_coli/SOM/parameters.txt", keep_local=True),
    output:
        "irsom/model/Escherichia_coli/SLSOM/W.txt",
        "irsom/model/Escherichia_coli/SLSOM/biases.txt",
        "irsom/model/Escherichia_coli/SLSOM/parameters.txt",
        "irsom/model/Escherichia_coli/SOM/units.txt",
        "irsom/model/Escherichia_coli/SOM/parameters.txt"
    shell:
        """
        mkdir -p irsom/model/Escherichia_coli/
        mv forge.ibisc.univ-evry.fr/lplaton/IRSOM/raw/master/model/Escherichia_coli/SLSOM irsom/model/Escherichia_coli/
        mv forge.ibisc.univ-evry.fr/lplaton/IRSOM/raw/master/model/Escherichia_coli/SOM irsom/model/Escherichia_coli/
        """

rule irsomGetFeaturer:
    input:
        HTTP.remote("https://forge.ibisc.univ-evry.fr/lplaton/IRSOM/raw/master/bin/Featurer", keep_local=True)
    output:
        "irsom/Featurer"
    shell:
        "mkdir -p irsom; mv {input} irsom/Featurer; chmod +x irsom/Featurer;"

rule predictIrsom:
   input:
       featurer="irsom/Featurer",
       slsomw="irsom/model/Escherichia_coli/SLSOM/W.txt",
       slsomb="irsom/model/Escherichia_coli/SLSOM/biases.txt",
       slsomp="irsom/model/Escherichia_coli/SLSOM/parameters.txt",
       somu="irsom/model/Escherichia_coli/SOM/units.txt",
       somp="irsom/model/Escherichia_coli/SOM/parameters.txt",
       transcripts="transcripts/{condition}-{replicate}/transcripts.fa",
   output:
       "irsom/{condition}-{replicate}/result.txt"
   threads: 1
   singularity:
        "docker://gelhausr/irsom:latest"
   shell:
       """
       mkdir -p irsom;
       predict.py --featurer={input.featurer} --file={input.transcripts} --model=irsom/model/Escherichia_coli/ --output=irsom/{wildcards.condition}-{wildcards.replicate}/
       """

rule irsomGFF:
    input:
        "irsom/{condition}-{replicate}/result.txt"
    output:
        "irsom/{condition, [a-zA-Z]+}-{replicate,\d+}.irsom.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/irsomGFF.py -c {wildcards.condition}  -i {input} -o {output}"

rule concatIrsom:
    input:
        lambda wildcards: expand("irsom/{{condition}}-{replicate}.irsom.gff", zip, replicate=samples.loc[(samples["method"] == "RNA") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.irsom.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; RiboReport/scripts/concatGFF.py {input} -o {output}"