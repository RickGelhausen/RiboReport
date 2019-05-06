rule pruneTranscripts:
    input:
        annotation="annotation/annotation.bed",
        genome="genomes/genome.fa"
    output:
        "riborater/transcripts.bed"
    conda:
        "../envs/riborater.yaml"
    threads: 20
    shell:
        "mkdir -p riborater; tools/ORF-RATER/prune_transcripts.py --inbed {input.annotation} --summarytable riborater/tid_removal_summary.txt -v {input.genome} --outbed {output} -p {threads}"

rule makeTfams:
    input:
        "riborater/transcripts.bed"
    output:
        "riborater/tfams/OUTSTEM.bed"
    conda:
        "../envs/riborater.yaml"
    threads: 1
    shell:
        "mkdir -p riborater; tools/ORF-RATER/make_tfams.py --inbed {input} --tfamstem riborater/tfam -v"

rule findOrfsAndTypes:
    input:
        annotation="riborater/transcripts.bed",
        tfam="riborater/tfams/OUTSTEM.bed",
        genome="genomes/genome.fa"
    output:
        "riborater/orfs.h5"
    conda:
        "../envs/riborater.yaml"
    threads: 20
    shell:
        "mkdir -p riborater; tools/ORF-RATER/find_orfs_and_types.py {input.genome} --inbed={input.annotation} --tfamstem riborater/tfam -p {threads} --codons ATG --orfstore {output} -v"

rule psiteTrimmed:
    input:
        annotation="annotation/annotation.bed",
        bam="maplink/RIBO-{condition}-{replicate}.bam"
    output:
        offsets="riborater/{condition}-{replicate}/offsets.txt",
        tally="riborater/{condition}-{replicate}/tallies.txt"
    conda:
        "../envs/riborater.yaml"
    threads: 20
    shell:
        "mkdir -p riborater; tools/ORF-RATER/psite_trimmed.py {input.bam} --minrdlen 15 --maxrdlen 29 --tallyfile {output.tally} --cdsbed {input.annotation} -p {threads} --offsetfile {output.offsets} -v"

rule regressOrfs:
    input:
        annotation="riborater/transcripts.bed",
        offsets="riborater/{condition}-{replicate}/offsets.txt",
        orfs="riborater/orfs.h5",
        bam="maplink/RIBO-{condition}-{replicate}.bam"
    output:
        regress="riborater/{condition}-{replicate}/regression.h5",
        metagene="riborater/{condition}-{replicate}/metagene.txt"
    conda:
        "../envs/riborater.yaml"
    threads: 20
    shell:
        "mkdir -p riborater; tools/ORF-RATER/regress_orfs.py {input.bam} --inbed {input.annotation} --offsetfile {input.offsets} --orfstore {input.orfs} --mincdsreads 10  -p {threads} --regressfile {output.regress} --metagenefile {output.metagene} -v"

rule rateRegressionOutput:
    input:
        regress="riborater/{condition}-{replicate}/regression.h5",
        orfs="riborater/orfs.h5"
    output:
        rating="riborater/{condition}-{replicate}/orfratings.h5",
        csv="riborater/{condition}-{replicate}/orfratings.csv"
    conda:
        "../envs/riborater.yaml"
    threads: 20
    shell:
        "mkdir -p riborater; tools/ORF-RATER/rate_regression_output.py --regressfile {input.regress} --orfstore {input.orfs} --CSV {output.csv} --ratingsfile {output.rating} -p {threads} -v"

rule makeOrfBed:
    input:
        annotation="riborater/transcripts.bed",
        ratings="riborater/{condition}-{replicate}/orfratings.h5"
    output:
        "riborater/{condition}-{replicate}/ratedorfs.bed"
    conda:
        "../envs/riborater.yaml"
    threads: 1
    shell:
        "mkdir -p riborater; tools/ORF-RATER/make_orf_bed.py --inbed {input.annotation} --ratingsfile {input.ratings} --outbed {output} --minrating 0"

rule quantifyOrfs:
    input:
        annotation="riborater/transcripts.bed",
        ratings="riborater/{condition}-{replicate}/orfratings.h5",
        offsets="riborater/{condition}-{replicate}/offsets.txt",
        metagene="riborater/{condition}-{replicate}/metagene.txt",
        bam="maplink/RIBO-{condition}-{replicate}.bam"
    output:
        qaunt="riborater/{condition}-{replicate}/quant.h5",
        csv="riborater/{condition}-{replicate}/quant.csv"
    conda:
        "../envs/riborater.yaml"
    threads: 20
    shell:
        "mkdir -p riborater; tools/ORF-RATER/quantify_orfs.py {input.bam} --inbed {input.annotation} --offsetfile {input.offsets} --metagene {input.metagene} --ratingsfile {input.ratings} --quantfile {output.quant} --CSV {output.csv} -p {threads} -v --minrating 0"
