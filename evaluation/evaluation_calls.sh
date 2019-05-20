#!/usr/bin/env bash

trap ctrl_c INT

function ctrl_c() {
echo "** Trapped CTRL-C"
exit
}

# Translatome vs tools
python3 statistics.py -r=data/all_filtered_massspec.gtf -t=data/combined.gtf -c=0.01 -o ./translatome_cov=0.01/ 
python3 statistics.py -r=data/all_filtered_massspec.gtf -t=data/combined.gtf -c=0.5 -o ./translatome_cov=0.5/ 

python3 plot_barplots.py -i1=./translatome_cov=0.01/df_stat.csv -i2=./translatome_cov=0.5/df_stat.csv -o=./translatome/ 

python3 venn_diagram.py -i=./translatome_cov=0.5/df_venn_genes.csv -o=./translatome/ -n=translatome -c=05
python3 venn_diagram.py -i=./translatome_cov=0.01/df_venn_genes.csv -o=./translatome/ -n=translatome -c=001

# SmallORFs: Translatome vs tools
python3 statistics.py -r=data/smallORFs.gtf -t=data/combined.gtf -c=0.01 -o ./smallORFs_cov=0.01/ 
python3 statistics.py -r=data/smallORFs.gtf -t=data/combined.gtf -c=0.5 -o ./smallORFs_cov=0.5/ 

python3 plot_barplots.py -i1=./smallORFs_cov=0.01/df_stat.csv -i2=./smallORFs_cov=0.5/df_stat.csv -o=./smallORFs/ 

python3 venn_diagram.py -i=./smallORFs_cov=0.5/df_venn_genes.csv -o=./smallORFs/ -n=smallORFs -c=05
python3 venn_diagram.py -i=./smallORFs_cov=0.01/df_venn_genes.csv -o=./smallORFs/ -n=smallORFs -c=001

# operon: Translatome vs tools
python3 statistics.py -r=data/ref_operon_intersect.gtf -t=data/tools_operon_intersect.gtf -c=0.01 -o ./operon_intersect_cov=0.01/ 
python3 statistics.py -r=data/ref_operon_intersect.gtf -t=data/tools_operon_intersect.gtf -c=0.5 -o ./operon_intersect_cov=0.5/ 

python3 plot_barplots.py -i1=./operon_intersect_cov=0.01/df_stat.csv -i2=./operon_intersect_cov=0.5/df_stat.csv -o=./operon_intersect/ 

python3 venn_diagram.py -i=./operon_intersect_cov=0.5/df_venn_genes.csv -o=./operon_intersect/ -n=operon_intersect -c=05
python3 venn_diagram.py -i=./operon_intersect_cov=0.01/df_venn_genes.csv -o=./operon_intersect/ -n=operon_intersect -c=001

# operon non overlapping: Translatome vs tools
python3 statistics.py -r=data/ref_operon_complement.gtf -t=data/tools_operon_complement.gtf -c=0.01 -o ./operon_complement_cov=0.01/ 
python3 statistics.py -r=data/ref_operon_complement.gtf -t=data/tools_operon_complement.gtf -c=0.5 -o ./operon_complement_cov=0.5/ 

python3 plot_barplots.py -i1=./operon_complement_cov=0.01/df_stat.csv -i2=./operon_complement_cov=0.5/df_stat.csv -o=./operon_complement/ 

python3 venn_diagram.py -i=./operon_complement_cov=0.5/df_venn_genes.csv -o=./operon_complement/ -n=operon_complement -c=05
python3 venn_diagram.py -i=./operon_complement_cov=0.01/df_venn_genes.csv -o=./operon_complement/ -n=operon_complement -c=001

# ncRNAs: Translatome vs tools
python3 statistics.py -r=data/ncRNAs.gtf -t=data/combined.gtf -c=0.01 -o ./ncRNAs_cov=0.01/ 
python3 statistics.py -r=data/ncRNAs.gtf -t=data/combined.gtf -c=0.5 -o ./ncRNAs_cov=0.5/ 

python3 plot_barplots.py -i1=./ncRNAs_cov=0.01/df_stat.csv -i2=./ncRNAs_cov=0.5/df_stat.csv -o=./ncRNAs/ 

python3 venn_diagram.py -i=./ncRNAs_cov=0.5/df_venn_genes.csv -o=./ncRNAs/ -n=ncRNAs -c=05
python3 venn_diagram.py -i=./ncRNAs_cov=0.01/df_venn_genes.csv -o=./ncRNAs/ -n=ncRNAs -c=001

# Pseudogenes: Translatome vs tools
python3 statistics.py -r=data/pseudogenes.gtf -t=data/combined.gtf -c=0.01 -o ./pseudogenes_cov=0.01/ 
python3 statistics.py -r=data/pseudogenes.gtf -t=data/combined.gtf -c=0.5 -o ./pseudogenes_cov=0.5/ 

python3 plot_barplots.py -i1=./pseudogenes_cov=0.01/df_stat.csv -i2=./pseudogenes_cov=0.5/df_stat.csv -o=./pseudogenes/ 

python3 venn_diagram.py -i=./pseudogenes_cov=0.5/df_venn_genes.csv -o=./pseudogenes/ -n=pseudogenes -c=05
python3 venn_diagram.py -i=./pseudogenes_cov=0.01/df_venn_genes.csv -o=./pseudogenes/ -n=pseudogenes -c=001

