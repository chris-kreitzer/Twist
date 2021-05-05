#!/bin/bash

files=( /scratch/kreitzer/tmp/*filtered_bcf )
length="${#files[@]}"

for ((RUN = 0; RUN < $length; RUN++));
do
        export RUN
        sbatch --export=ALL FilterVariants.sh
        sleep 1
done

