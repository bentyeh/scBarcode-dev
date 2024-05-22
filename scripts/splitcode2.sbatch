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
PATH_R2="$6"

splitcode -c "$PATH_CONFIG" \
	--nFastqs=2 --assign --mod-names --x-names --no-x-out --no-outb --out-fasta -t 4 \
	--mapping="$PATH_MAPPING" --summary="$PATH_SUMMARY" --output=/dev/stdout,/dev/null \
	"$PATH_R1" "$PATH_R2" |
awk '$1 ~ /^>/' |
pigz -p 4 > "$PATH_OUTPUT"