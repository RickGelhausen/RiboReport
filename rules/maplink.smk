rule maplink:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method}-{condition}-{replicate}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink; ln -s {params.inlink} {params.outlink}"

rule ribomaplink:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method, RIBO}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RIBO/; ln -s {params.inlink} {params.outlink}"

rule ribobamindexlink:
    input:
        "maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        "maplink/{method, RIBO}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.bai"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RIBO/; ln -s {params.inlink} {params.outlink}"

rule rnamaplink:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method, RNA}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RNA/; ln -s {params.inlink} {params.outlink}"

rule rnabamindexlink:
    input:
        "bam/{method}-{condition}-{replicate}.bam.bai"
    output:
        "maplink/{method, RNA}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.bai"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RNA/; ln -s {params.inlink} {params.outlink}"
