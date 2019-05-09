rule pruneTranscripts:
    input:
        annotation="annotation/annotation.bed",
        genome="genomes/genome.fa",
        bam="maplink/RIBO-{condition}-{replicate}.bam"
    output:
        "orfrater/{condition}-{replicate}/transcripts.bed"
    conda:
        "../envs/orfrater.yaml"
    threads: 20
    shell:
        "mkdir -p orfrater; ./tools/ORF-RATER/prune_transcripts.py  {input.genome} {input.bam} --inbed {input.annotation} --summarytable orfrater/{wildcards.condition}-{wildcards.replicate}/tid_removal_summary.txt -v --outbed {output} -p {threads}"

rule makeTfams:
    input:
        "orfrater/{condition}-{replicate}/transcripts.bed"
    output:
        "orfrater/{condition}-{replicate}/tfams.bed"
    conda:
        "../envs/orfrater.yaml"
    threads: 1
    shell:
        "mkdir -p orfrater; python tools/ORF-RATER/make_tfams.py --inbed {input} --tfamstem orfrater/{wildcards.condition}-{wildcards.replicate}/tfams -v"

rule findOrfsAndTypes:
    input:
        annotation="orfrater/{condition}-{replicate}/transcripts.bed",
        tfam="orfrater/{condition}-{replicate}/tfams.bed",
        genome="genomes/genome.fa"
    output:
        "orfrater/{condition}-{replicate}/orfs.h5"
    conda:
        "../envs/orfrater.yaml"
    threads: 20
    shell:
        "mkdir -p orfrater; python tools/ORF-RATER/find_orfs_and_types.py {input.genome} --inbed={input.annotation} --tfamstem orfrater/{wildcards.condition}-{wildcards.replicate}/tfams -p {threads} --codons ATG --orfstore {output} -v"

rule psiteTrimmed:
    input:
        #annotation="annotation/annotation.bed",
        annotation="orfrater/{condition}-{replicate}/transcripts.bed",
        bam="maplink/RIBO-{condition}-{replicate}.bam"
    output:
        offsets="orfrater/{condition}-{replicate}/offsets.txt",
        tally="orfrater/{condition}-{replicate}/tallies.txt"
    conda:
        "../envs/orfrater.yaml"
    threads: 20
    shell:
        "mkdir -p orfrater; python tools/ORF-RATER/psite_trimmed.py {input.bam} --minrdlen 22 --maxrdlen 40 --tallyfile {output.tally} --cdsbed {input.annotation} -p {threads} --offsetfile {output.offsets} -v"

rule regressOrfs:
    input:
        annotation="orfrater/{condition}-{replicate}/transcripts.bed",
        offsets="orfrater/{condition}-{replicate}/offsets.txt",
        orfs="orfrater/{condition}-{replicate}/orfs.h5",
        bam="maplink/RIBO-{condition}-{replicate}.bam"
    output:
        regress="orfrater/{condition}-{replicate}/regression.h5",
        metagene="orfrater/{condition}-{replicate}/metagene.txt"
    conda:
        "../envs/orfrater.yaml"
    threads: 20
    shell:
        "mkdir -p orfrater; python tools/ORF-RATER/regress_orfs.py {input.bam} --inbed {input.annotation} --offsetfile {input.offsets} --orfstore {input.orfs} --mincdsreads 0  -p {threads} --regressfile {output.regress} --metagenefile {output.metagene} -v"

rule rateRegressionOutput:
    input:
        regress="orfrater/{condition}-{replicate}/regression.h5",
        orfs="orfrater/{condition}-{replicate}/orfs.h5"
    output:
        rating="orfrater/{condition}-{replicate}/orfratings.h5",
        csv="orfrater/{condition}-{replicate}/orfratings.csv"
    conda:
        "../envs/orfrater.yaml"
    threads: 20
    shell:
        "mkdir -p orfrater; python tools/ORF-RATER/rate_regression_output.py {input.regress} --orfstore {input.orfs} --CSV {output.csv} --ratingsfile {output.rating} -p {threads} -v"

rule makeOrfBed:
    input:
        annotation="orfrater/{condition}-{replicate}/transcripts.bed",
        ratings="orfrater/{condition}-{replicate}/orfratings.h5"
    output:
        "orfrater/{condition}-{replicate}/ratedorfs.bed"
    conda:
        "../envs/orfrater.yaml"
    threads: 1
    shell:
        "mkdir -p orfrater; python tools/ORF-RATER/make_orf_bed.py --inbed {input.annotation} --ratingsfile {input.ratings} --outbed {output} --minrating 0"

rule quantifyOrfs:
    input:
        annotation="orfrater/{condition}-{replicate}/transcripts.bed",
        ratings="orfrater/{condition}-{replicate}/orfratings.h5",
        offsets="orfrater/{condition}-{replicate}/offsets.txt",
        metagene="orfrater/{condition}-{replicate}/metagene.txt",
        bam="maplink/RIBO-{condition}-{replicate}.bam"
    output:
        qaunt="orfrater/{condition}-{replicate}/quant.h5",
        csv="orfrater/{condition}-{replicate}/quant.csv"
    conda:
        "../envs/orfrater.yaml"
    threads: 20
    shell:
        "mkdir -p orfrater; python tools/ORF-RATER/quantify_orfs.py {input.bam} --inbed {input.annotation} --offsetfile {input.offsets} --metagene {input.metagene} --ratingsfile {input.ratings} --quantfile {output.quant} --CSV {output.csv} -p {threads} -v --minrating 0"
