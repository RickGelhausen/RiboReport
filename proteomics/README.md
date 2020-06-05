# Proteomics

Mass spectrometry data was retrieved from the supplemental materials of published data and converted to .gff format to allow comparison with the benchmark set.

Following are the exact conversion steps for each dataset:

## Escherichia coli:

Published proteomics data [Schmidt et al](https://doi.org/10.1038/nbt.3418) was obtained from [Supplementary Table S9](https://static-content.springer.com/esm/art\%3A10.1038\%2Fnbt.3418/MediaObjects/41587_2016_BFnbt3418_MOESM18_ESM.xlsx) by filtering column C with non-zero values of column Q of the cited manuscript.

Dataset2 was selected, because it is based on multiple biological replicates (see Figure S4)

```
Table S5 saved as tsv Nature_biotech_e_coli_proteome_tables.tsv

Extract Gene names and Proteins/cell column for LB 

  cut -f 2,7 Nature_biotech_e_coli_proteome_tables.tsv | tail -n +2 > ecoli_massspec_all.tsv
  
Remove empty lines and values that have 0 proteins/cell remove second column afterwards

  cat ecoli_massspec_all.tsv | grep -v -P "\t0" > ecoli_massspec_nonzero.tsv
  
  cut -f1 ecoli_massspec_nonzero.tsv > ecoli_massspec.tsv

  ./filter.py --in_gff_filepath ../data/escherichia_coli/annotation.gtf --feature_type "CDS" --filter_tag gene --in_locus_tag_filepath escherichia_coli/ecoli_massspec.tsv > ecoli_masspec.gff
```

## Pseudomonas aeroginosa:

Published proteomics data [Grady et al](https://dx.doi.org/10.1186%2Fs12864-017-3708-4) was obtained from Supplementary table [S21-S24](https://dx.doi.org/10.1186%2Fs12864-017-3708-4) of the cited manuscript.

```
Table saved as pseudomonas_aeroginosa/pseudomonas_aeroginosa_masspec.tsv

  ./filter.py --filter_tag "locus_tag" --feature_type "CDS" --in_gff_filepath ../data/pseudomonas_aeroginosa/annotation.gtf --in_locus_tag_filepath pseudomonas_aeroginosa/pseudomonas_aeroginosa_masspec.tsv > pseudomonas_aeroginosa_masspec.gff
```

## Listeria monocytogenes:

Published proteomics data [Impens et al](https://doi.org/10.1038/nmicrobiol.2017.5) was obtained from Supplementary Tables
[S2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S2.xlsx),
[S3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S3.xlsx),
[S4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S4.xlsx),
[S5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S5.xlsx),
[S6](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S6.xlsx),
[S7](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S7.xlsx),
[S8](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S8.xlsx)

```
# Retrieve massspec tables
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S2.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S3.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S4.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S5.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S6.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S7.xlsx
  
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802382/bin/NIHMS75486-supplement-Supplementary_table_S8.xlsx
  
# Convert all to csv

  libreoffice --headless --convert-to csv --outdir somedir *.xls

# Get name colum, drop header

  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S2.csv | tail -n +3 > S2.tsv
  
  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S3.csv | tail -n +3 > S3.tsv
  
  cut -d ',' -f2 NIHMS75486-supplement-Supplementary_table_S4.csv | tail -n +3 > S4.tsv
  
  cut -d ',' -f2 NIHMS75486-supplement-Supplementary_table_S5.csv | tail -n +3 > S5.tsv
  
  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S6.csv | tail -n +3 > S6.tsv
  
  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S7.csv | tail -n +3 > S7.tsv
  
  cut -d ',' -f1 NIHMS75486-supplement-Supplementary_table_S8.csv | tail -n +3 > S8.tsv
  
# S4 and S5 needed removal of inline comments, empty lines, S2 needs removal of empty line
#concat all tsv

  cat S*.tsv >> listeria_monocytogenes_massspec.tsv

# Table saved as: 

  ./filter.py --feature_type "gene"  --filter_tag "locus_tag" --in_gff_filepath ../data/listeria_monocytogenes/annotation.gff --in_locus_tag_filepath listeria_monocytogenes/listeria_monocytogenes_massspec.tsv > listeria_monocytogenes_masspec.gff
```

## Salmonella enterica serovar Typhimurium strain 14028s:

Published proteomics data [Yoon et al](https://dx.doi.org/10.1186%2F1752-0509-5-100) was obtained from Supplementary table [S1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3213010/bin/1752-0509-5-100-S1.XLSX) of the cited manuscript.

```
# Extracted first column with protein ids currently used by NCBI

# Changed identifiers to be compatible

  sed -i 's/STM/STM14_/g' 1752-0509-5-100-S1.tsv

# Removed entries only representing peptides
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

  ./masspec_from_locus_tags.py -a ../data/salmonella_enterica/annotation.gtf -m salmonella_enterica/1752-0509-5-100-S1.tsv -o salmonella_enterica/salmonella_masspec.gff
```
