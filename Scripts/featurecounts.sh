#!/bin/bash

featureCounts -T 16 -t exon -g gene_id -a /scratch/jmontenegro/nvectensis/results/annotation/tcs.gtf -o /proj/ferrer/rna_seq_twist/X204SC20120808-Z01-F001/raw_data/results/counts
*.bam

done


