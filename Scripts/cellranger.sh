#!/bin/bash

#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --output=/scratch/kreitzer/Nematostella/cellranger_%j.log
#SBATCH --error=/scratch/kreitzer/Nematostella/cellranger_%j.err
#SBATCH --export=ALL

# ENVIRONMENT
module load cellranger
module list

# CONSTANTS
wd="/scratch/kreitzer/Nematostella"
reads="/proj/agcole/data/rawdata_scRNAseq/Nematostella/tissue.pharynx.twistcontrol"
sample="87675"
res="${wd}/results"
od="${res}/map/cellranger"
ref="/proj/agcole/NvDb/Nv.cupcake.genome.new"

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

echo "cellranger count --id=pharynx.twistcontrol --fastqs=${sample} --sample=${sample} --transcriptome=Nv.cupcake.genome.new --nosecondary --force-cells=50000 --localcores=20"
cellranger count --id=pharynx.twistcontrol --fastqs=${sample} --sample=${sample} --transcriptome=Nv.cupcake.genome.new --nosecondary --force-cells=50000 --localcores=20

# move results to final directory
echo "mkdir ${od}"
mkdir ${od}

echo "cp -r pharynx.twistcontrol ${od}"
cp -r pharynx.twistcontrol ${od}

echo "Finished at `date`"
