#!/bin/bash

#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --output=/scratch/jmontenegro/alison/aaurita/logs/cellranger_%j.log
#SBATCH --error=/scratch/jmontenegro/alison/aaurita/logs/cellranger_%j.err
#SBATCH --export=ALL

# ENVIRONMENT
module load cellranger
module list

# CONSTANTS
wd="/scratch/jmontenegro/alison/aaurita"
reads="/proj/agcole/agcole/rawdata/Aurelia/96477"
sample="96477"
res="${wd}/results"
od="${res}/map/cellranger"
ref="${wd}/results/combine/aaurita.dedup"

# VARIABLES
### none

# EXECUTION
echo "Started at `date`"

# copy files to local disk
echo "cp -r ${reads} ${TMPDIR}"
cp -r ${reads} ${TMPDIR}
echo "cp -r ${ref} ${TMPDIR}"
cp -r ${ref} ${TMPDIR}

echo "cd ${TMPDIR}"
cd ${TMPDIR}

echo "cellranger count --id=polyp_TR2_50000 --fastqs=${sample} --sample=${sample} --transcriptome=aaurita.dedup --nosecondary --force-cells=50000 --localcores=20"
cellranger count --id=polyp_TR2_50000 --fastqs=${sample} --sample=${sample} --transcriptome=aaurita.dedup --nosecondary --force-cells=50000 --localcores=20

# move results to final directory
echo "mkdir -p ${od}"
mkdir -p ${od}

echo "cp -r polyp_TR2_50000 ${od}"
cp -r polyp_TR2_50000 ${od}

echo "Finished at `date`"
