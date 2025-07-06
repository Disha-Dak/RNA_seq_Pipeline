#!/bin/bash/

SECONDS = 0

cd /Users/dishajain/Desktop/RNA_seq_analysis
#Download sequence 
esearch -db sra -query PRJNA258216 \
  | efetch -format runinfo \
  > runinfo.csv

#Get wget of seq 
cut -d',' -f10 runinfo.csv \
  | tail -n +2 \
  | xargs -n1 wget -P raw_fastq/
#Step 1:Quality control- Fastqc
fastqc raw_fastq

java -jar /Users/dishajain/Desktop/RNA_seq_analysis SE -threads 4 raw_fastq_trimmed.fastq TRAILING:10 -phred33
echo "Trimmomatic finished running!"

fastqc raw_fastq_trimmed.fastq

# mkdir HISAT2
# get the genome indices
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

# run alignment
hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U raw_fastq_trimmed.fastq | samtools sort -o HISAT2/raw_fastq_trimmed.bam
echo "HISAT2 finished running!"


# get gtf
# wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
featureCounts -S 2 -a ../hg38/Homo_sapiens.GRCh38.106.gtf -o quants/demo_featurecounts.txt HISAT2/raw_fastq_trimmed.bam
echo "featureCounts finished running!"


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."


