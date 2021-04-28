#!/bin/bash

for RUN in $(seq 1 15);
do
	export RUN
	sbatch --export=ALL VariantCalling.sh
	sleep 1
done
	
