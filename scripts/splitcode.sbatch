#!/bin/bash
#
#SBATCH --time=01:00:00
#SBATCH -c 4
#SBATCH --mem-per-cpu=4G

source ~/.bashrc

PATH_CONFIG="$1"
PATH_MAPPING="$2"
PATH_SUMMARY="$3"
PATH_OUTPUT="$4"
PATH_R1="$5"

splitcode -c "$PATH_CONFIG" \
	--nFastqs=1 --assign --mod-names --x-names --no-x-out --no-outb --out-fasta -t 8 \
	--mapping="$PATH_MAPPING" --summary="$PATH_SUMMARY" --output=/dev/stdout \
	"$PATH_R1" |
awk '$1 ~ /^>/' |
pigz -p 8 > "$PATH_OUTPUT"