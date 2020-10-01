# Benchmark datasets

## Sequencing data sources:
Retived 4 different datasets for the construction of our benchmark dataset.
1. E. coli
[Ribo-seq and RNA-seq: GSE131514](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131514)
2. Listeria monocytogenes
[Ribo-seq: SAMEA3864955](https://www.ncbi.nlm.nih.gov/biosample/SAMEA3864955) and [RNA-seq: SAMEA3864956](https://www.ncbi.nlm.nih.gov/biosample/SAMEA3864956) datasets
3. Pseudomonas aeruginosa 
[Ribo-seq and RNA-seq: SAMN06617371](www.ncbi.nlm.nih.gov/biosample/SAMN06617371)
4. Salmonella enterica serovar Typhimurium
[Ribo-seq](https://www.ncbi.nlm.nih.gov/sra/SRX3456030[accn]) and [RNA-seq](https://www.ncbi.nlm.nih.gov/sra/SRX3456038[accn])

## Labeling
Using the annotation, all annotated ORFs were manually inspected and labeled as translated or not translated for each Ribo-seq dataset. The manual labels provide use with a ground truth dataset which can be used to benchmark Ribo-seq tools.


## Data for each organism
Each of the organims folder contains:
- annotation.gtf: NCBI annotation used for manual labeling
- config.yaml: the pipeline config file containing the adapter sequence
- genome.fa: genome file used for mapping
- labels_[organism].gff: all genes of the annotation but including the labeling information, expr=+ or expr=-,  in coloum 9 
- samples.tsv: the pipeline sample file containing the samples used for each organism
- misc_[organism].zip: contaning the preprocessed labeling information of the Ribo-Seq and the MS data. During the processing, different subsets were generated. 
These different sets are: CDS, ORFs part of operons (operons_intersect), ORFs not part of operons (operons_complement), sORFs (samllORFs)

### Details for the content of the misc_[organism].zip
- [subset]_labels_neg.gff
- [subset]_labels_pos.gff
- [subset]_masspec_pos.gff
- [organim]_masspec_pos.gff
- labels.gff
- labels_sorted.gff
- predicitons.gff
- smallORFs.txt
