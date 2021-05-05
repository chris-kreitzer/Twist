#!/bin/bash

files=( /proj/ferrer/rna_seq_twist/X204SC20120808-Z01-F001/raw_data/results/map/*bam )
length="${#files[@]}"

for ((RUN = 0; RUN < $length; RUN++));
do
	export RUN
	sbatch --export=ALL VariantCallingSample.sh
	sleep 1
done
