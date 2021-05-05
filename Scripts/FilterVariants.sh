#!/bin/bash

#SBATCH --job-name=Filtering
#SBATCH --cpus-per-task=6
#SBATCH --mem=4GB
#SBATCH --output=/scratch/kreitzer/logs/Filtering_%j.log
#SBATCH --error=/scratch/kreitzer/logs/Filtering_%j.err
#SBATCH --export=ALL

# use bcftools query for min.VAF = 0.2; and select specific position; twist variants

module load vcftools

# Environment
wd="/scratch/kreitzer"
tmp="/scratch/kreitzer/tmp/"

# Variables
echo "RUN = $RUN"
bcfiles=( ${tmp}/*filtered_bcf ) 
files="${bcfiles[$RUN]}"
base=`basename ${files}`
baseout="${base%.R1*}"


# Processing
cd ${tmp}

# filter variants 500 bp up- and downstream of twist gene
vcftools --bcf ${files} --out ${baseout}_twist_bcf --chr chr2 --from-bp 2357926 --to-bp 2360015 --recode --recode-INFO-all



