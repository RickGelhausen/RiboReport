# RiboReport19
This repository contains all code to recreate the contents of "RiboReport19 - Benchmarking ribosome profiling based identification of open reading frames in bacteria. Following are descriptions of the required steps.

## Processing of High Throughput Sequencing Data

## Generation of the Dataset

### Dependencies
- miniconda3
- snakemake =5.4.5

### Input data
For running the workflow, several input files are required:
- genome.fa (in the data folder)
- annotation.gtf (in the data folder)
- fastq files (Download link provided as soon as the GEO upload is finished)

### Workflow
---
**Note**
Even though snakemake workflows are executable locally, we do not advise this due to high memory usage and runtime of some of the processing steps. We ran the workflow on a TORQUE cluster system and provide according configuration files for TORQUE and SQE.

---


To run the provided snakemake workflow, follow the instructions below:

#### 1. Setup the workflow folder and download the workflow:

   ~~~~
   mkdir benchmark; cd benchmark;
   git clone git@github.com:RickGelhausen/ribo_benchmark.git
   ~~~~

#### 2. Fetch the annotation and genome files:

   ~~~~
   cp ribo_benchmark/data/annotation.zip . ;
   unzip annotation.zip;
   cp ribo_benchmark/data/genome.zip . ;
   unzip genome.zip;
   ~~~~

#### 3. Retrieve the sequencing data:
---
**Note**
This section will be updated as soon as the data GEO upload is finished.

---

#### 4. Run the snakemake workflow:

    In order to run snakemake, the creation of a conda environment is required. First install [miniconda3](https://docs.conda.io/en/latest/miniconda.html).

    Once miniconda3 is installed. Create a snakemake environment:
    ~~~~
    conda create -n snakemake -c conda-forge -c bioconda snakemake
    ~~~~

    Then you can copy and complete one of the provided submission scripts, or create your own.
    ~~~~
    cp ribo_benchmark/torque.sh .
    ~~~~
    or
    ~~~~
    cp ribo_benchmark/sge.sh .
    ~~~~

    Example for torque.sh:

    ~~~~
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
    snakemake --latency-wait 600 --use-conda -s ribo_benchmark/Snakefile --configfile ribo_benchmark/config.yaml --directory ${PWD} -j 20 --cluster-config ribo_benchmark/torque.yaml --cluster "qsub -N {cluster.jobname} -S /bin/bash -q {cluster.qname} -d <file path>/benchmark -l {cluster.resources} -o {cluster.logoutputdir} -j oe"
    ~~~~

    All **file path** statements have to be replaced by the path to your benchmark folder.

    **Please note** that these scripts might need some extra changes depending on your cluster system. If your cluster system does not run SGE or TORQUE, these files will most likely not work at all. In the case they run one SGE or TORQUE, there might be slightly different definitions for the resource statements (here `#PBS -l nodes=1:ppn=1`). This is then also the case for the configuration files `sge.yaml` and `torque.yaml`.


## Visualisation and plotting of the Results
All data can be visualised using the following scripts inside the evaluation folder. Further, an `evaluation_call.sh` script, containing all calls for the plotting pipeline, is provided. This bash script only works if all python scripts are positioned in the same folder and the input gtf-files from the data folder are stored in a "data"-folder in the same location.

### Dependencies
- numpy =1.16.3
- matplotlib =3.0.3
- seaborn =0.9.0
- pandas =0.24.2
- simple_venn =0.1.0 

### statistics.py

Parameters:
- reference_data (-r)
- tool_data (-t)
- save_path (-o)
- overlap_cutoff (-c)

Output:
- df_venn_FN_gene_dict.csv
- df_stat.csv
- df_venn_FP_predictions_dict.csv
- df_venn_predictions.csv
- df_venn_genes.csv


### plot_barplots.py
Parameters:
- input1_df (-i1)
- input2_df (-i2)
- save_path (-o)

Output:
- bar_FNR.pdf
- bar_recall.pdf
- bar_precision.pdf
- bar_FDR.pdf
- bar_F1.pdf

### venn_diagram.py
Parameters:
- input_df (-i)
- save_path (-o)
- name_folder (-n)
- coverage_percent (-c)

Output:
- venn_diagram.pdf 

