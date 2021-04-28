#!/bin/bash

#SBATCH --job-name=bulkRNA.variantCalls
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --output=/scratch/kreitzer/logs/bulkRNA.variantCalls_%j.log
#SBATCH --error=/scratch/kreitzer/logs/bulkRNA.variantCalls_%j.err
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

# Processing:
echo "Started at `date`"

<<<<<<< HEAD:Scripts/VariantCalling.sh
echo "bcftools mpileup -f ${ref} ${bam} | bcftools call -mv -Ob -o ${wd}_bcf"

bcftools mpileup -f ${ref} ${bam} | bcftools call -mv -Ob -o ${wd}_bcf

=======
for i in ${AlignedSamples[@]};
do
	name=$(basename $i)  \
	bcftools mpileup -f ${ref} $i | bcftools call -mv -Ob -o "/scratch/kreitzer/$name.bcf"  \
done
>>>>>>> 1c8706fcc9e20ccb58ce11485820cb0962c57557:Scripts/bulk.variantcall.sh

echo "Finished at `date`"

