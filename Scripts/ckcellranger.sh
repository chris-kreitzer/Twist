#!/bin/bash

#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --output=/scratch/kreitzer/Nematostella/cellranger2_%j.log
#SBATCH --error=/scratch/kreitzer/Nematostella/cellranger2_%j.err
#SBATCH --export=ALL

# ENVIRONMENT
module load cellranger
module list

# CONSTANTS
wd="/scratch/kreitzer/Nematostella"
ref="/proj/agcole/NvDb/Nv.cupcake.genome.new"
reads="/proj/agcole/data/rawdata_scRNAseq/Nematostella/"
fastqs="/proj/agcole/data/rawdata_scRNAseq/Nematostella/"
res="${wd}/results"
od="${res}/map/cellranger2"


# VARIABLES
echo "RUN = $RUN"
samples=(tissue.pharynx.twistcontrol tissue.pharynx.twistmutant)
readsprocess=${reads}${samples[$RUN]}
fastqsprocess=${fastqs}${samples[$RUN]}

cd ${readsprocess}
sample=`ls -h | head -n 1`
sampleprocess=${sample:0:5}
sampleID=`echo ${samples[$RUN]} | tr '.' '_'`

# EXECUTION
echo "Started at `date`"

# copy files to local disk
echo "cp -r ${readsprocess} ${TMPDIR}"
cp -r ${readsprocess} ${TMPDIR}
echo "cp -r ${ref} ${TMPDIR}"
cp -r ${ref} ${TMPDIR}

echo "cd ${TMPDIR}"
cd ${TMPDIR}

echo "cellranger count --id=${sampleID} --fastqs=${fastqsprocess} --sample=${sampleprocess} --transcriptome=Nv.cupcake.genome.new --nosecondary --force-cells=50000 --localcores=16"
cellranger count --id=${sampleID} --fastqs=${fastqsprocess} --sample=${sampleprocess} --transcriptome=Nv.cupcake.genome.new --nosecondary --force-cells=50000 --localcores=16

# move results to final directory
echo "mkdir ${od}"
mkdir ${od}

echo "cp -r ${sampleID} ${od}"
cp -r ${sampleID} ${od}

echo "Finished at `date`"
