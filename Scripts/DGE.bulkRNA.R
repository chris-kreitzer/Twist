## Differential gene expression from bulk RNA-Seq (Patricios experiment):
## I used featureCounts to get the reads for every feature; 
## now I will move on with DESeq2 to obtain a list of genes, which is differentially expressed in mutants and wild-type animals
## We need this list for further analysis in scRNA data (somehow as a sanity check)
## 
## Start: 05/20/2021
## chris kreitzer


setwd('~/Documents/GitHub/Twist/')
rm(list = ls())
.rs.restartR()


## Libraries

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2", force = T)
library("DESeq2")



## Input and Processing
bulkRNA_in = read.csv(file = 'Data_out/bulkRNAcounts_featureCounts', 
                      skip = 1, 
                      sep = '\t')
row.names(bulkRNA_in) = bulkRNA_in$Geneid
bulkRNA_in = bulkRNA_in[, seq(7, ncol(bulkRNA_in), 1)]

# set condition for DESeq2 analysis (sample basis)
sampleRNA = factor(c('Bubble1', 'Bubble2', 'Bubble3', 'Twi4d_1', 'Twi4d_2', 'Twi4d_3', 'Twihead_1', 'Twihead_2',
              'Twihead_3', 'WT4d_1', 'WT4d_2', 'WT4d_3', 'WThead_1', 'WThead_2', 'WThead_3'))

condition = factor(c('mut', 'mut', 'mut', 'mut', 'mut', 'mut', 
                     'mut', 'mut', 'mut', 'wt','wt', 'wt', 
                     'wt', 'wt', 'wt'))

# creating DESeq2 object:
bulkRNA_object = DESeqDataSetFromMatrix(countData = bulkRNA_in, 
                                    colData = DataFrame(sampleRNA, condition),
                                    design = sampleRNA ~ condition)


## DESeq2 WORKFLOW:
#' Pre-filtering the dataset:
#' remove rows with 0's: NO COUNTS OR just single counts across all samples
nrow(bulkRNA_object)
keep = rowSums(counts(bulkRNA_object)) > 1
bulkRNA_object = bulkRNA_object[keep, ]


#' The variance stabilizing transformation and the rlog
# Many common statistical methods for exploratory analysis of multidimensional data, 
# (clustering and PCA), work best for data that generally has the same range of variance at 
# different ranges of the mean values.

#' here I am using Variance-stabilizing transformation (VST);
vsd_bulkRNA = vst(bulkRNA_object, blind = FALSE)



head(assay(vsd), 3)



