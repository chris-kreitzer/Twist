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
bamfiles=( /scratch/kreitzer/data/*.bam )

# Variables
echo "RUN = $RUN"

chromosomes=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 )
chromosomeRUN="${chromosomes[$RUN]}"
base=`basename ${chromosomeRUN}`
baseout="${base%.R1*}"

# Processing:
echo "Started at `date`"

cd ${wd}
echo "bcftools mpileup -f ${ref} -r ${chromosomeRUN} -q 20 ${bamfiles[@]} | bcftools call -mv | bcftools +fill-tags -Oz -o tmp/${baseout}_bcf -- -t AF"
bcftools mpileup -f ${ref} -r ${chromosomeRUN} -q 20 ${bamfiles[@]} | bcftools call -mv | bcftools +fill-tags -Oz -- -t AF | bcftools filter -i 'AF>0.2' -Oz -o tmp/${baseout}_filtered_bcf
  

echo "Finished at `date`"


