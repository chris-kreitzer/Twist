#!/bin/bash

chromosomes=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15  )
length="${#chromosomes[@]}"

for ((RUN = 0; RUN < $length; RUN++));
do
  	export RUN
        sbatch --export=ALL VariantCallingChromosome.sh
        sleep 1
done



