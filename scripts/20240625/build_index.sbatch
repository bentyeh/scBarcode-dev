#!/bin/bash
#SBATCH -c 12
#SBATCH --mem-per-cpu=6000
#SBATCH --time=12:00:00

source ~/.bashrc
conda activate chipdip

PATH_HG38='/central/scratch/btyeh/hg38'
PATH_MM10='/central/scratch/btyeh/mm10'
PATH_INDEX='/central/scratch/btyeh/index_hg38_mm10'
mkdir -p "$PATH_INDEX"

bowtie2-build --threads 12 \
    -f "$(ls "$PATH_HG38"/h_chr*.fa "$PATH_MM10"/m_chr*.fa | paste -sd "," -)" \
    "$PATH_INDEX"/hg38_mm10
