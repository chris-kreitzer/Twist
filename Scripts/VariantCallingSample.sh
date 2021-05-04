#!/bin/bash

#SBATCH --job-name=VariantCallingRNA
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --output=/scratch/kreitzer/logs/VariantCalls_%j.log
#SBATCH --error=/scratch/kreitzer/logs/VariantCalls_%j.err
#SBATCH --export=ALL

# Setup:
module load bcftools

# Environment
wd="/scratch/kreitzer/"
ref="/scratch/jmontenegro/nvectensis/data/refs/nv_dovetail_4_gapped_chroms.final.fasta"
projfolder="/proj/ferrer/rna_seq_twist/X204SC20120808-Z01-F001/raw_data/results/map"

# Variables
echo "RUN = $RUN"
reads=( ${projfolder}/*.out.bam )
bam="${reads[$RUN]}"
base=`basename ${bam}`
baseout="${base%.R1*}"

# Processing:
echo "Started at `date`"

cd ${wd}
echo "bcftools mpileup -f ${ref} ${bam} | bcftools call -mv | bcftools +fill-tags -Oz -- -t AF | bcftools filter -i 'AF>0.2' -o tmp/${baseout}_filtered_bcf"

bcftools mpileup -f ${ref} ${bam} | bcftools call -mv | bcftools +fill-tags -Oz -- -t AF | bcftools filter -i 'AF>0.2' -Oz -o tmp/${baseout}_filtered_bcf

echo "Finished at `date`"

