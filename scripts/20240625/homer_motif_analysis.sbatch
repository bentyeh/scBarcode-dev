#!/bin/bash
#SBATCH -c 12
#SBATCH --mem-per-cpu=3000
#SBATCH --time=12:00:00

source ~/.bashrc
findMotifsGenome.pl "$PATH_BED" "$GENOME" "$DIR_ANALYSIS" -size "$SIZE" -p 12 &> "$PATH_LOG"
