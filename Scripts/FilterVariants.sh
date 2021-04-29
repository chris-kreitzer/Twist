#!/bin/bash

#SBATCH --job-name=Conversion
#SBATCH --cpus-per-task=6
#SBATCH --mem=4GB
#SBATCH --output=/scratch/kreitzer/logs/conversion_%j.log
#SBATCH --error=/scratch/kreitzer/logs/conversion_%j.err
#SBATCH --export=ALL

# use vcftools to filter for min.VAF = 0.2 (maf); and select specific position; twist variants

module load bcftools
module load vcftools

# Environment
wd="/scratch/kreitzer"
tmp="/scratch/kreitzer/tmp/"

# Variables
echo "RUN = $RUN"
bcfiles=( ${wd}/*_bcf ) 
files="${bcfiles[$RUN]}"
base=`basename ${files}`
baseout="${base%.R1*}"


# Processing
cd ${tmp}

# format the entire variants; --maf 0.2 (discuss)
vcftools --bcf ${files} --out ${baseout}_all --recode --recode-INFO-all

# filter variants 1000 bp up- and downstream of twist gene
vcftools --bcf ${files} --out ${baseout}_twist --chr chr2 --from-bp 2357626 --to-bp 2360015 --recode --recode-INFO-all




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
# first: bcf to vcf conversion
# bcftools view *_bcf > _vcf

# second: add AF (allele frequency info tag) to VCF output
# add only AF
# bcftools +fill-tags file.vcf  -- -t AF > output.vcf

# third: filter use bcftools query to filter variants with <= 0.2 AF
# bcftools query -i'INFO/AF>0.2' -f '%POS %REF %ALT %INFO\n' file.bcf
