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

## Labaling
Using the annotiaion, a ribo-seq dataset and its assosiated manually labled all ORFs as translated or not translated. Due to the mannual labels we 
obtained a ground truth dataset to which can be used for the benchmarking of Ribo-seq tools.


## Data for each organism
Each of the organims folder contains:
- annotaion.gtf: NCBI annotaion used for mannual labeling
- config.yaml
- genome.fa: genome file used for mapping
- labels_[organism].gff: all genes of the annotaion but including the lable information, expr=+ or expr=-,  in coloum 9 
- samples.tsv
- misc_[organism].zip: contaning the preprocessed lable information of the Ribo-Seq and the MS data. While the processing different subsets are generated. 
This diffrent sets are: CDS, ORFs part of operons (operons_intersect), ORFs not part of operons (operons_complement), sORFs (samllORFs)

### Deails content of the misc_[organism].zip
- [subset]_lables_neg.gff
- [subset]_lables_pos.gff
- [subset]_masspec_pos.gff
- [organim]_masspec_pos.gff
- labels.gff
- labels_sorted.gff
- predicitons.gff
- smallRFs.txt
