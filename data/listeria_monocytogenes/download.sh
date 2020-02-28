source activate sra-tools

fasterq-dump ERR1248499; pigz -p 2 ERR1248499.fastq; mv ERR1248499.fastq.gz LM-RNAseq-mRNA-Lm.fastq.gz
fasterq-dump ERR1248497; pigz -p 2 ERR1248497.fastq; mv ERR1248497.fastq.gz LM-RiboSeq-Lm.fastq.gz

mv *.fastq.gz fastq/

source deactivate
