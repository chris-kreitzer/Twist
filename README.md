# Twist - analysis of twist mutant in *Nematostella vectensis*

## Data logistics and setup:  

- **reference genome Nv:** /scratch/jmontenegro/nvectensis/data/refs/nv_dovetail_4_gapped_chroms.final.fasta.gz  
- **functional annotation NV:** /scratch/jmontenegro/nvectensis/data/orthologs/orth.table.tsv (gene predicted with corresponding functional annotation; alike BLAST query)   
- **twist RNA-seq data (bulk RNA Seq):**  
-  /proj/ferrer/rna_seq_twist   
-  /proj/ferrer/rna_seq_twist/X204SC20120808-Z01-F001/raw_data/results/map **[trimmed; aligned and sorted .bam files];** 
-  /proj/ferrer/rna_seq_twist/resequence/X204SC20120808-Z01-F004/raw_data/results/map **[re-sequenced samples];**  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Bubble (twist mutant with phenotype)  
TwiHead (twist mutant without phenotype)  
WTHead (as the name says)  
Twi4d (twist mutant 4dpf)  
WT4d (WT 4dpf)     
All were done as triplicates and have the suffix 1,2,3 added to their name depending on the replicate   

For the resequenced libraries the have the RS suffix added to them.   
Bubble 1,2 & 3   
Twihead 2  
WTHead 1 & 2  
Compared were:  
Bubble vs TwiHead  
Bubble vs WT  
Twihead vs WT  
Twi4d vs WT    
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Background:
**RNA-Seq** (bulk) were qc'ed, trimmed and aligned to the nemVec2 genome using STAR. STAR is currently the best INTRON-aware aligner;

**twist**; *Nv2.10864*; chr2: 2,358,626 - 2,359,015 (gene coordinates and location: 389 bp long); [mRNA](https://simrbase.stowers.org/feature/Nematostella/vectensis/mRNA/NVEC200_000161_1.1)

## First issues: 04/26/2021  
### bcftools [default seetings]:
bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf   

[/scratch/kreitzer]: bcftools mpileup -f /scratch/jmontenegro/nvectensis/data/refs/nv_dovetail_4_gapped_chroms.final.fasta  /proj/ferrer/rna_seq_twist/X204SC20120808-Z01-F001/raw_data/results/map/[...Aligned.bam] | bcftools call -mv -Ob -o calls.test.bcf  

one working example on /scratch/



## Generals:  
https://wiki.csb.univie.ac.at/doku.php?id=content:working_environment:services:slurm_workload_manager   
Working is intented to take place on /scratch
Processing and analysis should be conducted there. Once pipeline is through (/scratch/) ~ move data to **/proj/**
