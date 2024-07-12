#!/bin/bash
#SBATCH -c 10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=04:00:00

DIR_PROC="/central/groups/guttman/btyeh/scBarcode/data_proc/20240625"
DIR_TRIMMED="$DIR_PROC/trimmed"
BOWTIE2_PREFIX="/central/scratch/btyeh/index_hg38_mm10/hg38_mm10"

source ~/.bashrc
conda activate chipdip

bowtie2 \
    -p 10 \
    -t \
    --phred33 \
    -x "$BOWTIE2_PREFIX" \
    --maxins 2500 \
    -1 "$DIR_TRIMMED"/genomic_R1_val_1.fq.gz \
    -2 "$DIR_TRIMMED"/genomic_R2_val_2.fq.gz |
samtools view -b -o "$DIR_PROC/genomic.bam"
