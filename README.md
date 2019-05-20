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

Even though snakemake workflows are executable locally, we do not advise this due to high memory usage and runtime of some of the processing steps.
We ran the workflow on a TORQUE cluster system and provide according configuration files for TORQUE and SQE.
---

To run the provided snakemake workflow, follow the instructions below:

1. Setup the workflow folder and download the workflow:

   ~~~~
   mkdir benchmark; cd benchmark;
   git clone git@github.com:RickGelhausen/ribo_benchmark.git
   ~~~~

2. Fetch the annotation and genome files:

   ~~~~
   cp ribo_benchmark/data/annotation.zip . ;
   unzip annotation.zip;
   cp ribo_benchmark/data/genome.zip . ;
   unzip genome.zip;
   ~~~~

3. Retrieve the sequencing data:
---
**Note**
This section will be updated as soon as the data GEO upload is finished.

---

4. Run the snakemake workflow:

    In order to run snakemake, we suggest the creation of a conda environment.



## Visualisation and plotting of the Results
All data can be visualised using the following scripts inside the evaluation folder. Further, an evaluation_call.sh script, containing all calls for the plotting pipeline, is provided. This bash script only works if all python scripts are positioned in the same folder and the input gtf-files from the data folder are stored in a "data"-folder in the same location.

### Dependencies
- numpy =1.16.3
- matplotlib =3.0.3
- seaborn =0.9.0
- pandas =0.24.2
- simple_venn =0.1.0 

### statistics.py

Parameters:
- reference_data (-r) path to the gtf file soreing all genes of invertigated genome or the to be inverstigated subset of genes
- tool_data (-t) path to gtf file containing all predichtions for the tools: deepribom, ribotish, reparation and irsom
- save_path (-o) path to a directory, where the result tables should be stored
- overlap_cutoff (-c) dieserd sequence overlap cutoff between gene and prediction

Output:
- df_stat.csv main reault talbel soring all computed statistical measurments
- df_venn_FN_gene_dict.csv
- df_venn_FP_predictions_dict.csv
- df_venn_predictions.csv
- df_venn_genes.csv

This scipt generats for a reference and tool prediteion gtf file several statisitcal measuments. First true positives (TP), false positives (FP) and false negatives (FN) are predicted. One prediction will be assoiated with one gene and counted as one true positve. The assosiation selection is based on the lowest p-value (0.05). All genes fulfilling a overlap cutoff will be counted as suboptimals and not as false positvies. False positives are all predictions, where no gene fulfilling the overlap cutoff could be found and the false negatives vice verser. Based on this computations the recall, FNR, precision, FDR and F1 measure are calculated. 

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

