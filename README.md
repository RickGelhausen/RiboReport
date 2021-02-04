<img src="RiboReport.png" width="661">

# RiboReport21
This repository contains all code to recreate the contents of "RiboReport '20 - Benchmarking ribosome-profiling based identification of open reading frames in bacteria". Following are descriptions of the required steps.

| WARNING: This documentation will be updated as soon as the sequencing data for the *escherichia_coli* is made public! |
| --- |

## Processing of High Throughput Sequencing Data
In the publication, the data required for the generation of our result figures was created using the `snakemake` workflow described below. After the creation of the input data, the figures from the publication were generated by using the evaluation scripts provided in the evaluation folder.
The generation and usage of the final tables and figures is described below.

## Generation of the Dataset

### Dependencies
- miniconda3
- snakemake =5.4.5

### Input data
For running the workflow, several input files are required:
- genome.fa (in the data folder)
- annotation.gtf (in the data folder)
- samples.tsv (in the data folder)
- config.yaml (in the data folder)
- fastq files (and bigwig files) available via NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131514, token:ofkbusuqzbcvfip)

### Workflow
---
**Note**
Even though `snakemake` workflows are executable locally, we do not advise this due to high memory usage and runtime of some of the processing steps. We ran the workflow on a TORQUE cluster system and provide according configuration files for TORQUE and SGE.

---

To run the provided `snakemake` workflow, follow the example for `salmonella_enterica` below:

#### 1. Setup the workflow folder and download the workflow:

~~~~
mkdir salmonella_enterica; cd salmonella_enterica;
wget https://github.com/RickGelhausen/RiboReport/archive/2020.tar.gz;
tar -xzf 2020.tar.gz; mv RiboReport-2020 RiboReport; rm 2020.tar.gz;
~~~~

#### 2. Required tools
In order to run the workflow, the tools analyzed in the publication have to be installed.
Reparation and RiboTISH are automatically downloaded from bioconda and do not need any prior installation.
IRSOM and DeepRibo are not on conda and have to be installed manually.

First of all, create a `tools` folder in the `salmonella_enterica` repository.

~~~~
mkdir tools; cd tools;
~~~~

Next, we install `IRSOM`. To make it compatible with our workflow, we first create a conda environment for the dependencies:

~~~~
conda create -n irsom -c bioconda -c conda-forge plotnine pandas numpy tensorflow matplotlib docopt python=3.6.8
conda activate irsom
~~~~
(use source activate, if conda is not set-up for your bash)

Then, we download and install `IRSOM`:

~~~~
git clone https://forge.ibisc.univ-evry.fr/lplaton/IRSOM.git
pip install -r IRSOM/pip_package.txt
~~~~

This should only install `IRSOM`, as the other dependencies are already installed in the conda environment.

~~~~
conda activate irsom
~~~~

This leaves the installation of DeepRibo. For DeepRibo all dependencies will be downloaded automatically by the workflow.
We just have to download DeepRibo itself.

~~~~
wget https://github.com/Biobix/DeepRibo/archive/v1.1.zip
unzip v1.1.zip; mv DeepRibo-1.1 DeepRibo; rm v1.1.zip; mv DeepRibo tools;
~~~~

#### 3. Fetch the annotation and genome files:

~~~~
cp RiboReport/data/salmonella_enterica/annotation.gtf . ;
cp RiboReport/data/salmonella_enterica/genome.fa . ;
~~~~

#### 4. Retrieve the sequencing data:

There are many ways to download fastq files with SRA. For more information about downloading please have a look at the following guide: [Downloading SRA data using command line utilities](https://www.ncbi.nlm.nih.gov/books/NBK158899/).

The simplest way is most likely the usage of the [SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std), as it allows direct conversion into `.fastq` files.

Using the `SRA Toolkit` and the `SRR ID` of each sample, we can use the `fasterq-dump` executable to download the according `.fastq` files.

If you do not have the `SRA Toolkit`, we suggest using the conda environment:

~~~~
conda create -n sra-tools -c bioconda -c conda-forge sra-tools pigz
conda activate sra-tools
~~~~
(use source activate, if conda is not set-up for your bash)

For simplicity, we already collected the required `fasterq-dump` calls in a file.

~~~~
cp RiboReport/data/salmonella_enterica/download.sh .
bash download.sh
~~~~

This will download all required fastq files into a fastq folder.

#### 5. Run the snakemake workflow:

In order to run `snakemake`, the creation of a conda environment is required. First install [miniconda3](https://docs.conda.io/en/latest/miniconda.html).

Once miniconda3 is installed. Create a snakemake environment:
~~~~
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
~~~~

Then you can copy and complete one of the provided submission scripts, or create your own.
~~~~
cp RiboReport/torque.sh .
~~~~
or
~~~~
cp RiboReport/sge.sh .
~~~~

Example for `torque.sh`:

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
snakemake --latency-wait 600 --use-conda -s RiboReport/Snakefile --configfile RiboReport/config.yaml --directory ${PWD} -j 20 --cluster-config RiboReport/torque.yaml --cluster "qsub -N {cluster.jobname} -S /bin/bash -q {cluster.qname} -d <file path>/benchmark -l {cluster.resources} -o {cluster.logoutputdir} -j oe"
~~~~

All **file path** statements have to be replaced by the path to your benchmark folder.

**Please note** that these scripts might need some extra changes depending on your cluster system. If your cluster system does not run SGE or TORQUE, these files will most likely not work at all. In the case they run SGE or TORQUE, there might be slightly different definitions for the resource statements (here `#PBS -l nodes=1:ppn=1`). This is then also the case for the configuration files `sge.yaml` and `torque.yaml`.


## Visualisation and plotting of the Results
All data can be visualised using the following scripts inside the evaluation folder. Further, an `evaluation_call.sh` script, containing all calls for the plotting pipeline, is provided. This bash script only works if all python scripts are positioned in the same folder and the input gtf-files from the data folder are stored in a "data"-folder in the same location.

### Dependencies
- numpy =1.16.3
- matplotlib =3.0.3
- seaborn =0.9.0
- pandas =0.24.2
- simple_venn =0.1.0
- bedtools =v2.28.0

### statistics_pos_neg.py

Parameters:
- reference_pos_data (-p) path to the .gtf file containing all genes of the investigated genome or the investigated subset of genes
- reference_neg_data (-n) path to the .gtf file containing all genes of the investigated genome or the investigated subset of genes
- tool_data (-t) path to gtf file containing all predictions for the tools: deepribo, ribotish, reparation and irsom
- save_path (-o) path to a directory, where the result tables will be stored
- overlap_cutoff (-c) allowed sequence overlap cutoff between gene and prediction
- flag_subopt (-s) allowed allow suboptimals in statistical calculation. Set to 1 to enable flag. If not 0
- flag_no_gene (-g) add predictions which do not overlap with any genen for the given cutoff. Set to 1 to enable flag. If not 0.

Output:
- df_stat.csv -> main result table storing all computed statistical measurements. Needed for the bar plots.
- SetX_labels_pos_predictions_overlap_fp.gtf -> all prdictions which do not overlap with the positive or the negative gene set. Therefore, predictions outside of the annotation borders for a given overlap (temporary result).
- SetX_labels_pos_predictions_overlap_neg.gtf -> all genes of the negative set overlapping with a given percent with any prediction (temporary result).
- SetX_labels_pos_predictions_overlap_pos.gtf -> all genes of the positive set overlapping with a given percent with any prediction (temporary result).
- SetX_labels_pos_predictions_overlap_temp.gtf -> all prdictions which do not overlap with the positive gene set (temporary result).
- df_venn_genes.csv -> table containing a column for each tool listing the number of genes counted as true positive (TP). Needed for the overlap plots.
- toolX_score_list -> For each prediction, correct genen prediction and not predicted gene a tuple of: (score, overlap, label) is saved as a python object serialization (pickle.dump). Needed to compute the ROC and PRC.

This script generates several statistical measurements for a reference and tool prediction .gtf file. First true positives (TP), false positives (FP) and false negatives (FN) are predicted. One prediction will be associated with one gene and counted as one true positve. The association selection is based on the lowest p-value (0.05). All genes fulfilling the overlap cutoff will be counted as suboptimals and not as false positives. False positives are all predictions, where no gene fulfilling the overlap cutoff could be found and the false negatives vice versa. Based on this computations the recall, FNR, precision, FDR and F1 measure are calculated.

### plot_barplots.py
Parameters:
- input1_df (-i1) df_stat.csv of the statistics.py script for the first overlap cutoff.
- input2_df (-i2) df_stat.csv of the statistics.py script for the second overlap cutoff.
- save_path (-o) path to directory where the plots will be stored.

Output:
- bar_FNR.pdf -> barplot of FNR measures for different tools and two cutoff conditions.
- bar_recall.pdf -> barplot of recall measures for different tools and two cutoff conditions.
- bar_precision.pdf -> barplot of precision measures for different tools and two cutoff conditions.
- bar_FDR.pdf -> barplot of FDR measures for different tools and two cutoff conditions.
- bar_F1.pdf -> barplot of F1 measures for different tools and two cutoff conditions.

This script generates barplots of statistical measures for the different tools and two overlap cutoff conditions. The output barplots consist of FNR, recall, precision, FDR and F1 measures. It is based on the statistics.py statistical output table.

### venn_diagram.py
Parameters:
- input_df (-i) .cvs containing the list of true positive genes for the used tools.
- save_path (-o) path where the Venn diagramm will be saved.
- name_folder (-n) name of the result folder of the investigated set. It will be used in the file name to make it unique.
- coverage_percent (-c) percentage of overlap cutoff that was used to determine the true positves. It will be used in the file name to make it unique.

Output:
- venn_diagram.pdf -> a 4 Venn diagram highlighting overlapping predictions of the tool.

This script will generate a 4 Venn diagram for the overlap of true positive predicted genes of the 4 investigated tools.

### prc_plotting.py
Parameters:
- input_path: path list of triplets (score, overlap, label) for each tool.
- species_label: name of the bacterial datasets.
- save_path: path to save data to.
- percent_overlap: how much of the prediction should overlap with the gene and vice versa.
- experiment: name of the subset.

Output:
- datesetX_prc_overlapScore_subset_lables_.pdf

This scripts generates Precision-Recall-Curves (PRC) for all different tools and overlaps.

### roc_plotting.py
Parameters:
- input_path: path list of triplets (score, overlap, label) for each tool.
- save_path: path to save data to.
- percent_overlap: how much of the prediction should overlap with the gene and vice versa.
- experiment: name of the subset.

Output:
- datesetX_roc_overlapScore_.pdf

This scripts generates Receiver-Operator-Curves (ROC) for all different tools and overlaps.

