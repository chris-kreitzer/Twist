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
# BiocManager::install("apeglm")
library("DESeq2")
library(ggplot2)
library(dplyr)
library("pheatmap")
library("RColorBrewer")
library("ggbeeswarm")
library("apeglm")
library("genefilter")


## Input and Processing
bulkRNA_in = read.csv(file = 'Data_out/bulkRNAcounts_featureCounts', 
                      skip = 1, 
                      sep = '\t')
row.names(bulkRNA_in) = bulkRNA_in$Geneid
bulkRNA_in = bulkRNA_in[, seq(7, ncol(bulkRNA_in), 1)]
colnames(bulkRNA_in) = c('Bubble1', 'Bubble2', 'Bubble3', 'Twi4d_1', 'Twi4d_2', 'Twi4d_3', 'Twihead_1', 'Twihead_2', 'Twihead_3', 
                         'WT4d_1', 'WT4d_2', 'WT4d_3', 'WThead_1', 'WThead_2', 'WThead_3')


#' raw data summary for twist gene across all samples
twist.raw = bulkRNA_in[which(row.names(bulkRNA_in) == 'NV2.10864'), ]

all.out = data.frame()
for(i in seq(1, ncol(twist.raw), 3)){
  t.mean = apply(twist.raw[1, seq(i, i+2, 1)], 1, mean)
  t.sub = data.frame(group = colnames(twist.raw)[i],
                     mean.twist = t.mean)
  all.out = rbind(all.out, t.sub)
}

#' make plot
ggplot(all.out, aes(x = group, y = mean.twist)) + geom_bar(stat = 'identity') +
  labs(title = 'RAW counts; most basic MEAN summary over 3 replicates')



# set condition for DESeq2 analysis (sample basis)
sampleRNA = factor(c('Bubble1', 'Bubble2', 'Bubble3', 'Twi4d_1', 'Twi4d_2', 'Twi4d_3', 'Twihead_1', 'Twihead_2', 'Twihead_3', 
                     'WT4d_1', 'WT4d_2', 'WT4d_3', 'WThead_1', 'WThead_2', 'WThead_3'))

condition = factor(c('mut', 'mut', 'mut', 'mut', 'mut', 'mut', 
                     'mut', 'mut', 'mut', 'wt','wt', 'wt', 
                     'wt', 'wt', 'wt'))

# creating DESeq2 object:
bulkRNA_raw = DESeqDataSetFromMatrix(countData = bulkRNA_in, 
                                    colData = DataFrame(sampleRNA, condition),
                                    design = sampleRNA ~ condition)


## DESeq2 WORKFLOW:
#' Pre-filtering the dataset:
#' remove rows with 0's: NO COUNTS OR just single counts across all samples
nrow(bulkRNA_raw)
keep = rowSums(counts(bulkRNA_raw)) > 1
bulkRNA_object = bulkRNA_raw[keep, ]


#' The variance stabilizing transformation and the rlog
# Many common statistical methods for exploratory analysis of multidimensional data, 
# (clustering and PCA), work best for data that generally has the same range of variance at 
# different ranges of the mean values.

#' here I am using Variance-stabilizing transformation (VST);
vsd_bulkRNA = vst(bulkRNA_object, blind = FALSE)


# make a plot on the VST transformation:
modi.bulkRNA_object = estimateSizeFactors(bulkRNA_object)

transformed_data = bind_rows(
  as_data_frame(log2(counts(modi.bulkRNA_object, normalized = T)[, 1:15] + 1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd_bulkRNA)[, 1:15]) %>% 
    mutate(transformation = "vst"))

colnames(transformed_data)[1:2] = c("x", "y")  

lvls = c("log2(x + 1)", "vst")
transformed_data$transformation = factor(transformed_data$transformation, levels = lvls)

VST_bulkRNA.plot = ggplot(transformed_data, 
       aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

VST_bulkRNA.plot

#' variance stabilized transformation for counts data for Twist;
vsd.twist = assay(vsd_bulkRNA)
vsd.twist = data.frame(vsd.twist[which(row.names(vsd.twist) == 'NV2.10864'), ])
colnames(vsd.twist) = 'vsd'
vsd.twist$sample = row.names(vsd.twist)
vsd.twist$sample = substr(x = vsd.twist$sample, start = 1, stop = nchar(vsd.twist$sample) - 1)

#' make plot
ggplot(vsd.twist, aes(x = sample, y = vsd)) + geom_jitter(width = 0.1, size = 1.3) +
  theme_bw() +
  coord_fixed() +
  labs(x = '', y = 'variance-stabilization-transformation', title = 'vst-Twist expression')


#' Sample distances: how similar are samples 
sampleDist = dist(t(assay(vsd_bulkRNA)))

sampleDistMatrix = as.matrix(sampleDist)
rownames(sampleDistMatrix) = paste(vsd_bulkRNA$sampleRNA, vsd_bulkRNA$condition, sep = " - " )
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         col = colors)


## PCA: sample to sample distance:
plotPCA(vsd_bulkRNA, intgroup = c("condition"))
plotPCA(vsd_bulkRNA, intgroup = c("sampleRNA", 'condition'))

#' return the rawData
PCA_raw = plotPCA(vsd_bulkRNA, intgroup = c( "sampleRNA", "condition"), returnData = TRUE)
ggplot(PCA_raw, 
       aes(x = PC1, y = PC2, 
           color = sampleRNA, 
           shape = condition)) +
  geom_point(size =3) +
  coord_fixed() +
  scale_color_manual(values = c('Bubble1' = 'red',
                                'Bubble2' = 'red',
                                'Bubble3' = 'red',
                                'Twi4d_1' = 'blue',
                                'Twi4d_2' = 'blue',
                                'Twi4d_3' = 'blue',
                                'Twihead_1' = 'aquamarine3',
                                'Twihead_2' = 'aquamarine3',
                                'Twihead_3' = 'aquamarine3',
                                'WT4d_1' = 'brown',
                                'WT4d_2' = 'brown',
                                'WT4d_3' = 'brown',
                                'WThead_1' = 'orange',
                                'WThead_2' = 'orange',
                                'WThead_3' = 'orange')) +
  ggtitle("PCA with variance-stabilized transformed data")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Differential gene expression: 
DGE_bulkRNA = DESeq(bulkRNA_object)
DGE_out = results(DGE_bulkRNA)

#' raise the log2 fold change threshold, 
#' show more substantial changes due to treatment (CRISPR)
#' instead of lfcThreshold = 0
DGE_out_conservative = results(DGE_bulkRNA, lfcThreshold = 1)
table(DGE_out_conservative$padj < 0.1)

#' extracting genes with strongest lfc down-regulation signal
results_out = subset(DGE_out_conservative, padj < 0.1)
results_out = results_out[order(results_out$log2FoldChange), ]
results_out = as.data.frame(results_out, row.names = row.names(results_out))

#' example; of DE gene; NV2.618 (can this trend already be seen in vst data?)
vsd.618 = assay(vsd_bulkRNA)
vsd.618 = data.frame(vsd.618[which(row.names(vsd.618) == 'NV2.618'), ])
colnames(vsd.618) = 'vsd'
vsd.618$sample = row.names(vsd.618)
vsd.618$sample = substr(x = vsd.618$sample, start = 1, stop = nchar(vsd.618$sample) - 1)
vsd.618$group = c(rep('mut', 9), rep('wt', 6))

#' make plot
ggplot(vsd.618, aes(x = sample, y = vsd, color = group)) + 
  geom_jitter(width = 0.1, size = 1.3) +
  scale_color_manual(values = c('mut' = 'red',
                                'wt' = 'blue')) +
  theme_bw() +
  coord_fixed() +
  labs(x = '', y = 'variance-stabilization-transformation', title = 'vst-NV2.618 expression')

## make a visualization for the top gene (most different)
topDiffGene = rownames(DGE_out_conservative)[which.min(DGE_out_conservative$padj)]
plotCounts(DGE_bulkRNA, gene = topDiffGene, intgroup = 'condition')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Differential gene expression; precise look into the data
topGene = plotCounts(DGE_bulkRNA, gene = topDiffGene, intgroup = c('condition', 'sampleRNA'), returnData = T)

#' visualizaiton
topGene.plot = ggplot(topGene, aes(x = condition, y = count, color = sampleRNA)) +
  scale_y_log10() +  
  geom_beeswarm(cex = 3) +
  theme_bw() +
  labs(y = 'log10(counts)', title = paste0(topDiffGene))
topGene.plot


#' MA plot:
resultsNames(DGE_bulkRNA);
resMA = lfcShrink(DGE_bulkRNA, coef = 'condition_wt_vs_mut', type = 'apeglm')
plotMA(resMA)


#' gene clustering:
topVarGenes = head(order(rowVars(assay(vsd_bulkRNA)), decreasing = TRUE), 20)
mat = assay(vsd_bulkRNA)[topVarGenes, ]
mat = mat - rowMeans(mat)
anno = as.data.frame(colData(vsd_bulkRNA)[, c('sampleRNA', 'condition')])
pheatmap(mat, annotation_col = anno, main = 'top20 diff gene; [vst-counts]', fontsize = 12)

#' one example;
x = assay(vsd_bulkRNA)[row.names(assay(vsd_bulkRNA)) == 'NV2.5412', ]












