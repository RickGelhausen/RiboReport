rule ribotishAnnotation:
    input:
        annotation="auxiliary/featurecount/annotation.gtf",
        sizes="genomes/sizes.genome"
    output:
        "ribotish/annotation_processed.gtf"
    conda:
        "../envs/annotation.yaml"
    threads: 1
    shell:
        "mkdir -p ribotish; ribo_benchmark/scripts/createRiboTISHannotation.py -a {input.annotation}  --genome_sizes {input.sizes} --annotation_output {output}"
