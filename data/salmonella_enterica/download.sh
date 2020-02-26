source activate sra-tools

fasterq-dump SRR6359968; pigz -p 2 SRR6359968.fastq; mv SRR6359968.fastq.gz RIBO-WT-3.fastq.gz

fasterq-dump SRR6359976; pigz -p 2 SRR6359976.fastq; mv SRR6359976.fastq.gz RNA-WT-3.fastq.gz

mv *.fastq.gz fastq/
source deactivate
