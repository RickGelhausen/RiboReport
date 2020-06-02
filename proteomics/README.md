# Escherichia coli:

Dataset2 was selected, because it is based on multiple biological replicates (see Figure S4)

Table S5 saved as tsv Nature_biotech_e_coli_proteome_tables.tsv

Extract Gene names and Proteins/cell column for LB 

  cut -f 2,7 Nature_biotech_e_coli_proteome_tables.tsv | tail -n +2 > ecoli_massspec_all.tsv
  
Remove empty lines and values that have 0 proteins/cell remove second column afterwards

  cat ecoli_massspec_all.tsv | grep -v -P "\t0" > ecoli_massspec_nonzero.tsv
  
cut -f1 ecoli_massspec_nonzero.tsv > ecoli_massspec.tsv

rsync -auv fr_fe1017@login02.binac.uni-tuebingen.de:/beegfs/work/fr_fe1017/rt/SPP2002/bm01/annotation/annotation.gtf

./filter.py --in_gff_filepath escherichia_coli/annotation.gtf --feature_type "CDS" --filter_tag gene --in_locus_tag_filepath escherichia_coli/ecoli_massspec.tsv > ecoli_masspec.gff


# Pseudomonas aeroginosa:

Table saved as pseudomonas_aeroginosa/pseudomonas_aeroginosa_masspec.tsv
  rsync -auv fr_fe1017@login02.binac.uni-tuebingen.de:/beegfs/work/fr_fe1017/rt/SPP2002/bm06/annotation/annotation.gtf pseudomonas_aeroginosa/annotation.gtf

  ./filter.py --filter_tag "locus_tag" --feature_type "CDS" --in_gff_filepath pseudomonas_aeroginosa/annotation.gtf --in_locus_tag_filepath pseudomonas_aeroginosa/pseudomonas_aeroginosa_masspec.tsv > pseudomonas_aeroginosa_masspec.gff

# Listeria monocytogenes:

  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S2.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S3.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S4.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S5.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S6.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S7.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S8.xlsx
  
convert all to csv

  libreoffice --headless --convert-to csv --outdir somedir *.xls

get name colum, drop header

  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S2.csv | tail -n +3 > S2.tsv
  
  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S3.csv | tail -n +3 > S3.tsv
  
  cut -d ',' -f2 NIHMS75486-supplement-Supplementary_table_S4.csv | tail -n +3 > S4.tsv
  
  cut -d ',' -f2 NIHMS75486-supplement-Supplementary_table_S5.csv | tail -n +3 > S5.tsv
  
  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S6.csv | tail -n +3 > S6.tsv
  
  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S7.csv | tail -n +3 > S7.tsv
  
  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S8.csv | tail -n +3 > S8.tsv
  
#S4 and S5 needed removal of inline comments, empty lines, S2 needs removal of empty line
#concat all tsv

  cat S*.tsv >> listeria_monocytogenes_massspec.tsv

Table saved as: 

  rsync -auv fr_fe1017@login02.binac.uni-tuebingen.de:/beegfs/work/fr_fe1017/rt/SPP2002/bm03/annotation/annotation.gtf listeria_monocytogenes/annotation.gff

./filter.py --feature_type "gene"  --filter_tag "locus_tag" --in_gff_filepath listeria_monocytogenes/annotation.gff --in_locus_tag_filepath listeria_monocytogenes/listeria_monocytogenes_massspec.tsv > listeria_monocytogenes_masspec.gff


# Salmonella enterica serovar Typhimurium strain 14028s:

Extracted first column with present protein ids

changed identifiers to be compatible

sed -i 's/STM/STM14_/g' 1752-0509-5-100-S1.tsv

removed entries only representing peptides

PSLT011
PSLT028
PSLT031
PSLT037
PSLT038
PSLT039
PSLT046
PSLT052
PSLT053
PSLT103
