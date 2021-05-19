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

## First summary about the CRISPR-Cas9 Off-targets approach:   
From the bulk-RNA Seq we observed 4 variants which are mutant-specific. We exclusively concentrated on INDELs (because they are most likely to produce truncating proteins). We run the pipeline chromosome-wise, meaning that we get population estimates from the input samples (sum 15). After annotation and filtering, we have a list of 4 INDEL variants which are of particular interest (present in all 9 mutant bulk-RNA samples, and missing in the remaining 6 wildtype samples).   
- chr2:1059227	CCCCGAAC	CC	gene_id NV2.10722; transcript_id NV2.10722.1; exon_number 1; this particular variant has a **delta bp of 6**, meaning that 2 AS are not produced (missing), however we still do assume that a protein is produced (perhaps altered function).  
- chr2:2358816	GGTAACGT	GGT	gene_id NV2.10864; transcript_id NV2.10864.1; exon_number 1; this is actually **twist** mutation. We see a **5bp** deletion in all the mutants; since this is an odd number to 6, we assume that this causes a true ***frameshift mutation***, most likely that a pre-mature stop codon is introduced and hence the protein function is disabled completely.  
- chr2:3292703	CTTCCATTT	CTTT	gene_id NV2.10979; transcript_id NV2.10979.1; exon_number 2; Orthologe: ***LOC100177304 [Ciona intestinalis]***, RSPH1_MOUSE Radial spoke head 1 homolog; RSPH1 Radial spoke head 1 homolog; this particular gene is affected and MUST be considered in further analysis. As there is a **5 bps deletion** we equally assume that this variant causes a pre-mature stop-codon, and hence the protein function might be disabled.   
- chr8:15289229	CTGTGTGTGTGTGTGTGTGTGT	CTGTGTGTGTGTGTGT|CTGTGTGTGTGTGTGTGTGT	gene_id NV2.24244; transcript_id NV2.24244.1; exon_number 19; **ditto** as for first variant described in this list -- no particular effect apart from altered protein function;  

Moreover, I discovered **5 variants** which are unique to the wild-type samples (e.g. TwiHead, WT4d, etc.) and are missing in the CRISPR-mutants. Meaning that those variants acquired the mutation 'naturally' (*de novo*), and CRISPR doesn't effect those genes in the mutants. As we are not particularly interested in those, we are not concentrating on this oberservation for this project.   

From this point of view, we will consider NV2.10979 and obviously *twist* in the scRNA data. The questions we want to address are: 
- are those genes expressed in the scRNA data (sanity check);
- if they are expressed, in which cluster are they expressed;
- how is the expression changing - different from wild-types to mutants;   


Before, I move on to scRNA libraries (data), I will define a GOI - differentially expressed in the mutants and wild-types. This set of genes should serve as a sanity check, whether we see the same pattern in scRNA data. 


## Background:
**RNA-Seq** (bulk) were qc'ed, trimmed and aligned to the nemVec2 genome using STAR. STAR is currently the best INTRON-aware aligner;  
We are specifically interested in variants in the bulk-RNA data, as single-cell data can't be used to find off-targets. We are using the RNA-bulk data as a sanity check for the single-cell data. Biological material for both sequencing strategies (bulk, single-cell) comes from the same animals!  
We are running bcftools mpileup population-wise (meaning all 15 samples together) and loop through chromosomes. Because we are running bcftools population-wise we can calculate proper AF values. The amount of alternative alleles (within 15 samples) divided by reference (amount). AF < 0.2 are discared.  
We then move on to find intersects (**bedtools**); we concentrate just on tcs_exon_gtf - as we only want to see variants on exons (likely functional effect); biological effects of intronic variants are hard to predict!   

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
