#!/bin/bash
#$ -N benchmark
#$ -cwd 
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -o <file path>/benchmark
#$ -e <file path>/benchmark
#$ -j y
#$ -R y
export PATH="<file path>/miniconda3/bin/:$PATH"
cd <file path>/benchmark
source activate snakemake
snakemake --use-conda -p -s ribo_benchmark/Snakefile --configfile ribo_benchmark/config.yaml --directory ${PWD} -j 20 --cluster-config ribo_benchmark/sge.yaml --cluster "qsub -R y -N {cluster.jobname} -cwd -pe {cluster.parallelenvironment} -l {cluster.memory} -o {cluster.logoutputdir} -e {cluster.erroroutputdir} -j {cluster.joinlogs}" --latency-wait 60

