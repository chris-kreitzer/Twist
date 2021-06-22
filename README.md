# Twist - analysis of twist mutant in *Nematostella vectensis*   


## Hypothesis about ***NV-twist***   
The observed CRISPR/Cas9 ***twist*** knock-out phenotype (Bubble, failed tentacle formation) is novel in Nematostella vectensis. 

***NV-twist*** is involved in tentacle formation (broader context; muscle cell differentiation)   
***NV-twist*** as a master TF regulator is involved in initiating a cancer-like phenotype (first described in Nematostella vectensis)

## The main goal is hence:   
- see which genes are differentially expressed in ***twist*** mutants and wild type animals
- can we define a unique cell cluster, which elucidate the cell identity of this novel Bubble phenotype   
- can we define a gene-regulatory network in which ***twist*** is involved. We are specifically interested in signalling pathways and Transcription factors.   

## We need the money because:   
- We have already produced a lot of data (single-cell RNA expression) and urgently a data analysist (ideally two)   
- Based on those analysis (whether there is a link between cancer and twist, etc.) we would further go into validation studies   
- Parallel we would employ a post-doc specialising in Antibody production, since there is no antibody available for ***twist***. An AB would be the definitiv goal, as with this in hand, we can start a Chip-Seq experiment, where we can directly elucidate the binding domains of ***twist*** on the DNA and hence infer the direct targets of ***NV-twist***  
- We further need money, because we aim to publish in a high-ranking journal (if we can define cancer-like phenotype) and hence our reputation as lab as well as the contribution to a broader science community is given. 


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


## scRNA analysis with Seurat:    
#### Background:   
We have used cellranger to map the RNA reads to the reference genome and obtained a count-matrix per cell. Individual cells are captures in the cellranger output, through ***unique*** Barcode sequences (UMI's). This means, Seurat uses a count-matrix, where columns represent unique DNA-Sequences (UMI's = cells) and rows are features (in the provided annotation table we have roughly 20,000 features for Nematostella). For every cell and feature we have a count (produced by cellranger; and will be loaded in Seurat.   
#### Workflow Seurat:   
First we look into features (per cell), RNA counts and percentage MT. These metrics are used as a quality measure. Cells where the the feature count is below 250 are discarded (not informative) and cells with RNA-counts >20,000 are equally discarded - as this may be indicative of a duplicate cell sequenced (two cells in 10X rather than a single Transcriptome). We filter out those cells, apply a normalization and scaling approach `scale_factor = 5,000` [Alison's suggestion]. 
We then merge the two libraries (Mutant and Control) and define a list of highly variable genes (by default, Seurat uses `vst` AND `nfeatures = 2000`).  
We then run a PCA `runPCA` (first dimision reduction). Every dot basically represents a cell with it's transcriptome. This means, with `runPCA` we take the 2,000 variable features (among cells) and cluster cells (dimension reduced) with similar transcriptomic profiles [runPCA example](https://github.com/chris-kreitzer/Twist/blob/main/Figures/Pharynx_PCAs_2000features.pdf).   
What we can clearly see from the plot:   
- the overall transcriptome (single cells) from both the mutants and the the control are very much comparable   
- There are some spatial structures (in which the arrangement of clusters doesn't play an important role; just that several **distinct** clusters exist).   
- What's really interesting is, that one particular cell cluster stands out in the mutant library (pharynx_mutant);   
- now it's time to take a closer look into this distinct cell cluster; to see where those cells (with its transcriptome belong to, or which identity those cells show)   
#### Look into a specific feature (***NvTwist***) among the two libraries, and see spatial expression (normalized) in the data / cell-identities.   
[quantitative (normalized) NvTwist expression among cell clusters](https://github.com/chris-kreitzer/Twist/blob/main/Figures/NvTwist_expression%20among%20tissues.pdf). You can see that NvTwist is expressed in both libraries (at a varying degree). However, consider the different density. It seems like that NvTwist is more expressed in controls and in mutant libraries (especially, the left corner).   

#### Cluster cells - important parameter: resolution: The higher the resolution, the more clusters will be formed (and vice versa).   
`Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset`.   
[How to select the resolution](https://github.com/chris-kreitzer/Twist/blob/main/Figures/Clustering_Resolution_clusttree.pdf). You see we start off with 14 clusters; Looking into the top-row, we see that at **low resolution** we have direct trajectories. However, there is one cell-identity lineage which does not change at any resolution, meaning that **9** stays very much the same and show distinct features. Whereas cluster 0 gets a little bit blurry and splits in many smaller clusters (at higher resolution). Looking at low-level plots (getting sense of data); a low resolution is good to go. However, we are interested in getting the most discinct clusters - just looking @this [plot](https://github.com/chris-kreitzer/Twist/blob/main/Figures/Clustering_Resolution_clusttree.pdf), we are choosing `resolution = 0.8`, because there we are keeping the most sensitiviy and structure. 

#### Visualizing the Cluster's) - UMAP (based on clustering): BuildClusterTree: Phylogenetic Analysis of Identity Classes:   
When we look into the reduced and clustered ID's [Figure](https://github.com/chris-kreitzer/Twist/blob/main/Figures/Clusters_UMAP_reduced.pdf) we can clearly see that there are two features which distiguish the mutants from the controls. Firstly, we have a cell cluster **20**(turquois) which is apparent in the mutants and missing in the the controls.   
Secondly, we can also see that cell cluster **13** may be of interest, as this cluster is enriched for **NvTwist** in the control and largely missing in the mutants. At least, **NvTwist** is lower expressed in the mutant library at this cluster. We can also see that the overall clustering pattern looks 'clean' in the sense, that there is little dispersion of different cell-identities with other clusters. Apart from **23 AND 24**. Look also at this [Figure](https://github.com/chris-kreitzer/Twist/blob/main/Figures/CellType%20Distribution%20among%20Libraries.pdf). There you can clearly see that cluster **20** seems to exclusively be just available at the mutant library where it is missing in the control.  



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
