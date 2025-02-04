#!/bin/bash

PATH_A="$1"
PATH_B="$2"
shift 2
OPTIONS="$@"

source ~/.bashrc
conda activate chipdip

bedtools intersect -c "$@" -a "$PATH_A" -b "$PATH_B" |
cut -f 7 |
awk '{s += $1} END {print s}'
