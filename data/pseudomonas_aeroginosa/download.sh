source activate sra-tools

fasterq-dump SRR5356902; pigz -p 2 SRR5356902.fastq; mv SRR5356902.fastq.gz RNA-PAO1-nalk-1.fastq.gz
fasterq-dump SRR5356901; pigz -p 2 SRR5356901.fastq; mv SRR5356901.fastq.gz RIBO-PAO1-nalk-1.fastq.gz

mv *.fastq.gz fastq/

source deactivate sra-tools
