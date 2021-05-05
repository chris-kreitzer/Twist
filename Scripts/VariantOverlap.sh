#!/bin/bash

#SBATCH --job-name=VariantOverlap
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --output=/scratch/kreitzer/logs/Overlaps_%j.log
#SBATCH --error=/scratch/kreitzer/logs/Overlaps_%j.err
#SBATCH --export=ALL

# Setup:
module load bedtools

# Environment
tmp="/scratch/kreitzer/tmp"
ref="/scratch/jmontenegro/nvectensis/results/annotation/tcs_exons_only.gtf"

# Variables
echo "RUN = $RUN"
files=( ${tmp}/*filtered_bcf )
filesprocess="${files[$RUN]}"
base=`basename ${filesprocess}`
baseout="${base%.R1*}"

# Processing:
echo "Started at `date`"
cd ${tmp}

echo "bedtools intersect -f 0.6 -wa -wb -sorted -a ${filesprocess} -b ${ref} > tmp/${baseout}_intersected.txt" 
bedtools intersect -u -wa -wb -sorted -a ${filesprocess} -b ${ref} > ${baseout}_intersected.txt

echo "Finished at `date`"
