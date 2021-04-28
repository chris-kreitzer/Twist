#!/bin/bash

#SBATCH --job-name=FilterVariants
#SBATCH --cpus-per-task=6
#SBATCH --mem=6GB
#SBATCH --output=/scratch/kreitzer/logs/FilterVariants_%j.log
#SBATCH --error=/scratch/kreitzer/logs/FilterVariants_%j.err
#SBATCH --export=ALL

module load bcftools


# first: bcf to vcf conversion
bcftools view *_bcf > _vcf

# second: add AF (allele frequency info tag) to VCF output
# add only AF
bcftools +fill-tags file.vcf  -- -t AF > output.vcf

# third: filter use bcftools query to filter variants with <= 0.2 AF

