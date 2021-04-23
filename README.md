# Twist - analysis of twist mutant in *Nematostella vectensis*

## Data logistics and setup:  

- **reference genome Nv:** /scratch/jmontenegro/nvectensis/data/refs/nv_dovetail_4_gapped_chroms.final.fasta.gz  
- **functional annotation NV:** /scratch/jmontenegro/nvectensis/data/orthologs/orth.table.tsv (gene predicted with corresponding functional annotation; alike BLAST query)   
- **twist RNA-seq data (bulk RNA Seq):** /proj/ferrer/rna_seq_twist | /proj/ferrer/rna_seq_twist/X204SC20120808-Z01-F001/raw_data/results/map (first time) | /proj/ferrer/rna_seq_twist/resequence/X204SC20120808-Z01-F004/raw_data/results/map (re-sequenced samples)   
 
- sample annotation for .bam files:   

Bubble (twist mutant with phenotype)  
TwiHead (twist mutant without phenotype)  
WTHead (as the name says)  
Twi4d (twist mutant 4dpf)  
WT4d (WT 4dpf)     
All were done as triplicates and have the suffix 1,2,3 added to their name depending on the replicate   

For the resequenced libraries the have the RS suffix added to them.  
The resequenced libraries are:   
Bubble 1,2 & 3   
Twihead 2  
WTHead 1 & 2  
Compared were:  
Bubble vs TwiHead  
Bubble vs WT  
Twihead vs WT  
Twi4d vs WT    

## First issues: 


## Generals:
Working is intented to take place on /scratch
Processing and analysis should be conducted there. Once pipeline is through (/scratch/) ~ move data to **/proj/**
