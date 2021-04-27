#!/bin/bash

#SBATCH --job-name=bulkRNA.variantCalls
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --output=/scratch/kreitzer/logs/bulkRNA.variantCalls_%j.log
#SBATCH --error=/scratch/kreitzer/logs/bulkRNA.variantCalls_%j.err
#SBATCH --export=ALL

# Setup:
module load bcftools


wd="/scratch/kreitzer/"
projfolder="/proj/ferrer/rna_seq_twist/X204SC20120808-Z01-F001/raw_data/results/map/"
AlignedSamples=$(find ${projfolder} -type f -iname "Aligned.sortedByCoord.out.bam")

ref="/scratch/jmontenegro/nvectensis/data/refs/nv_dovetail_4_gapped_chroms.final.fasta"

# Execution:
echo "Started at `date`"

for i in $AlignedSamples;
do
	name=$(basename $i)
	bcftools mpileup -f ${ref} $i | bcftools call -mv -Ob -o echo "/scratch/kreitzer/$name.bcf"

echo "Finished at `date`"

## running with default parameters
## out: 04/27/2021
