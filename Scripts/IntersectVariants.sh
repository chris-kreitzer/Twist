#!/bin/bash

module load bedtools

bedtools intersect -a /scratch/jmontenegro/nvectensis/results/annotation/tcs.gtf - b example.bed -wa -c

