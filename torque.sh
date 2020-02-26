#!/bin/bash
#PBS -N benchmark
#PBS -S /bin/bash
#PBS -q "long"
#PBS -d <file path>/benchmark
#PBS -l nodes=1:ppn=1
#PBS -o <file path>/benchmark
#PBS -j oe
cd <file path>/benchmark
export PATH="<file path>/miniconda3/bin/:$PATH"
source activate snakemake
snakemake --latency-wait 600 --use-conda -s RiboReport/Snakefile --configfile RiboReport/config.yaml --directory ${PWD} -j 20 --cluster-config RiboReport/torque.yaml --cluster "qsub -N {cluster.jobname} -S /bin/bash -q {cluster.qname} -d <file path>/benchmark -l {cluster.resources} -o {cluster.logoutputdir} -j oe"
