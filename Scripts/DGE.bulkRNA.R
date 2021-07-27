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
library(rlist)
library(xlsx)
library(cowplot)
library(forcats)

## Input and Processing
load('NV2_annotation')
bulkRNA_in = read.csv(file = 'Data_out/bulkRNAcounts_featureCounts', 
                      skip = 1, 
                      sep = '\t')
row.names(bulkRNA_in) = bulkRNA_in$Geneid
bulkRNA_in = bulkRNA_in[, seq(7, ncol(bulkRNA_in), 1)]
colnames(bulkRNA_in) = c('Bubble1', 'Bubble2', 'Bubble3', 'Twi4d_1', 'Twi4d_2', 'Twi4d_3', 'Twihead_1', 'Twihead_2', 'Twihead_3', 
                         'WT4d_1', 'WT4d_2', 'WT4d_3', 'WThead_1', 'WThead_2', 'WThead_3')


#' gene filtering; @least: 10 reads per replicate - in two groups to keep gene
#' Loop through Libraries (at least 10 reads in every replicate), and then at least in two libraries to compare DGE
df = list()
for(i in seq(1, 15, by = 3)){
  print(i)
  data.subset = bulkRNA_in[, c(i, i+1, i+2)]
  data.subset$var = ifelse(apply(data.subset, 1, min) >= 10, 'keep', 'discard')
  data.subset$gene = row.names(data.subset)
  colnames(data.subset)[grep(pattern = 'var', colnames(data.subset))] = colnames(bulkRNA_in)[i]
  data.keep = data.subset[, c(ncol(data.subset) - 1, ncol(data.subset))]
  
  df[[i]] = data.keep
  rm(data.keep)
  rm(data.subset)
  
}

#' sub modies;
df = df[lengths(df) != 0] # exclude zero elements
df = list.cbind(df)
df[, 'gene' == names(df)] = NULL
df$final = ifelse(apply(df, 1, function(x) sum(x == 'keep') >= 1), 'keep', 'discard') #' library merge

#' which genes to keep
gene.ii = row.names(bulkRNA_in)[which(df$final == 'keep')]

#' updated raw counts
bulkRNA_modi = bulkRNA_in[which(row.names(bulkRNA_in) %in% gene.ii),, drop = F]

#' convert to matrix
bulkRNA_modi = as.matrix(bulkRNA_modi)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set condition for DESeq2 analysis (sample basis)
# As we are solely interested in within group comparisions, we don't need to specify
# any specific desing formula. As the name itself already contain the conditions (WT vs Mutant and Time)
# 
# I don't include individuals in the formula, as all the individuals show roughly the same expression;
coldata = DataFrame(phenotype = factor(rep(c('Bubble', 'Twi4d', 'TwiHead', 'WT4d', 'WTHead'), each = 3)))

# creating DESeq2 object:
bulkRNA_object_phenotype = DESeqDataSetFromMatrix(countData = bulkRNA_modi, 
                                    colData = coldata,
                                    design = ~phenotype)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## DESeq2 WORKFLOW:
#' working with bulkRNA_object

#' The variance stabilizing transformation and the rlog
# Many common statistical methods for exploratory analysis of multidimensional data, 
# (clustering and PCA), work best for data that generally has the same range of variance at 
# different ranges of the mean values.

#' here I am using Variance-stabilizing transformation (VST);
vsd_bulkRNA = vst(bulkRNA_object_phenotype, blind = FALSE)
rlog_bulkRNA = rlog(bulkRNA_object_phenotype, blind = FALSE)

#' compare the normalization methods on twist;
bulkRNA_object_phenotype = estimateSizeFactors(bulkRNA_object_phenotype)

transformed_data = bind_rows(
  as_data_frame(log2(counts(bulkRNA_object_phenotype, normalized = T)[, 1:15] + 1)) %>%
    mutate(transformation = "log2(x + 1)",
           gene = row.names(assay(bulkRNA_object))),
  as_data_frame(assay(vsd_bulkRNA)[, 1:15]) %>% 
    mutate(transformation = "vst",
           gene = row.names(assay(vsd_bulkRNA))),
  as_data_frame(assay(rlog_bulkRNA)[, 1:15]) %>% 
    mutate(transformation = "rlog",
           gene = row.names(assay(rlog_bulkRNA))))


#' extracting normalized counts for twist (3 independent methods): 
twist_transformed = transformed_data[which(transformed_data$gene == 'NV2.10864'), ]

twist.converted = data.frame()
for(i in 1:nrow(twist_transformed)){
  out = as.data.frame(t(twist_transformed[i, ]))
  out$method = out$V1[which(row.names(out) == 'transformation')]
  out$library = row.names(out)
  out$library = gsub('\\_[0-9]$', '', out$library)
  out = out[1:15, ]
  twist.converted = rbind(twist.converted, out)
}


#' make plot; comparing normalization methods:
twist.converted$library = ifelse(grepl(pattern = 'Bubble.*', twist.converted$library), 
                                 gsub('[0-9]$', '', twist.converted$library), twist.converted$library)

twist_bulkRNA.plot = ggplot(twist.converted,
       aes(x = library, y = as.numeric(V1))) + 
  geom_jitter(width = 0.2) +
  facet_grid( . ~ method) +
  theme_bw() +
  labs(x = '', y = 'normalized counts',
       title = 'Comparing different normalization methods on twist expression across different libraries')

twist_bulkRNA.plot

#' the main take-home from here is that different normalization methods do not skew the overall pattern;



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
plotPCA(vsd_bulkRNA, intgroup = c("phenotype"))

#' return the rawData
PCA_raw = plotPCA(vsd_bulkRNA, intgroup = c( "phenotype"), returnData = TRUE)
ggplot(PCA_raw, 
       aes(x = PC1, y = PC2, 
           color = phenotype)) +
  geom_point(size = 3) +
  coord_fixed() +
  scale_color_manual(values = c('Bubble' = 'red',
                                'Twi4d' = 'blue',
                                'TwiHead' = 'aquamarine3',
                                'WT4d' = 'brown',
                                'WTHead' = 'orange')) + 
  ggtitle("PCA on variance-stabilized transformed data") +
  labs(x = 'PC1: 82% variance', y = 'PC2: 12% variance') +
  guides(color = guide_legend(title = 'Pheno-/Genotype\ncombination')) +
  theme_bw() +
  theme(axis.text = element_blank())



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Differential gene expression: 
#' DESeq function consits of the following steps:
#' estimation of size factors (controlling for differences in the sequencing depth of the samples);
#' estimation of dispersion values for each gene
#' fitting a generalized linear model (negative binomial model)
DGE_bulkRNA = DESeq(bulkRNA_object_phenotype)

#' get differential expression results;
#' Bubble vs WT4d
res = results(object = DGE_bulkRNA, lfcThreshold = 1, contrast = c('phenotype', 'TwiHead', 'WTHead'))

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

#' make annotations for the DE genes within the groups; focus on TF and muscle development;
ortho_table = read.csv(file = 'Data_out/NV2_orthologe.table.tsv', sep = '\t')
ortho_table = ortho_table[!duplicated(ortho_table$NV2Id) & !is.na(ortho_table$NV2Id), ]

#' one example: Twist 4d (mutant) early development vs Bubble (adult animal with mutation and phenotype)
TwiHead_WTHead = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
          by.x = 'Gene', by.y = 'NV2Id', all.x = T)


#' refine based on p-adjust and logFC
TwiHead_WTHead = TwiHead_WTHead[which(abs(TwiHead_WTHead$log2FoldChange) > 1 & TwiHead_WTHead$padj < 0.01),, drop = F]
TwiHead_WTHead = TwiHead_WTHead[, -grep('^Bubble*', colnames(TwiHead_WTHead))]
TwiHead_WTHead = TwiHead_WTHead[, -grep('WT4d*', colnames(TwiHead_WTHead))]
TwiHead_WTHead = TwiHead_WTHead[, -grep('Twi4d*', colnames(TwiHead_WTHead))]

#' add meta data to up-/down table
TwiHead_WTHead_full = merge(x = TwiHead_WTHead, y = NV2_annotation, 
                            by.x = 'Gene', by.y = 'NV2', all.x = T)
write.xlsx(x = TwiHead_WTHead_full, file = 'Data_out/TwiHead_WTHead_fulltable.xlsx', row.names = F)
# TwiHead_WTHead_full$gene_short_name[which(TwiHead_WTHead_full$log2FoldChange < 0 & !is.na(TwiHead_WTHead_full$TF))]



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' comparision of Twi4d and WT4d; whether we see that twist is signigicantly different expressed;
res = results(object = DGE_bulkRNA, lfcThreshold = 1, contrast = c('phenotype', 'Twi4d', 'WT4d'))

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

Twi4d_WT4d = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
                    by.x = 'Gene', by.y = 'NV2Id', all.x = T)

#' refine based on padjust and logFC
Twi4d_WT4d = Twi4d_WT4d[which(abs(Twi4d_WT4d$log2FoldChange) > 1 & Twi4d_WT4d$padj < 0.01),, drop = F]
Twi4d_WT4d = Twi4d_WT4d[, -grep('^WThead*', colnames(Twi4d_WT4d))]
Twi4d_WT4d = Twi4d_WT4d[, -grep('Twihead*', colnames(Twi4d_WT4d))]
Twi4d_WT4d = Twi4d_WT4d[, -grep('Bubble*', colnames(Twi4d_WT4d))]

#' add meta data to up-/down table
Twi4d_WT4d_full = merge(x = Twi4d_WT4d, y = NV2_annotation, 
                            by.x = 'Gene', by.y = 'NV2', all.x = T)
write.xlsx(x = Twi4d_WT4d_full, file = 'Data_out/Twi4d_WT4d_fulltable.xlsx', row.names = F)
#' Twi4d_WT4d_full$gene_short_name[which(Twi4d_WT4d_full$log2FoldChange > 0 & !is.na(Twi4d_WT4d_full$TF))]


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' comparision of Bubble and Twihead; whether we see that twist is signigicantly different expressed;
res = results(object = DGE_bulkRNA, lfcThreshold = 1, contrast = c('phenotype', 'Bubble', 'TwiHead'))

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

Bubble_TwiHead = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
                   by.x = 'Gene', by.y = 'NV2Id', all.x = T)

#' refine based on padjust and logFC
Bubble_TwiHead = Bubble_TwiHead[which(abs(Bubble_TwiHead$log2FoldChange) > 1 & Bubble_TwiHead$padj < 0.01),, drop = F]
Bubble_TwiHead = Bubble_TwiHead[, -grep('^WThead*', colnames(Bubble_TwiHead))]
Bubble_TwiHead = Bubble_TwiHead[, -grep('Twi4d*', colnames(Bubble_TwiHead))]
Bubble_TwiHead = Bubble_TwiHead[, -grep('WT4d*', colnames(Bubble_TwiHead))]

#' add meta data to up-/down table
Bubble_TwiHead_full = merge(x = Bubble_TwiHead, y = NV2_annotation, 
                        by.x = 'Gene', by.y = 'NV2', all.x = T)
write.xlsx(x = Bubble_TwiHead_full, file = 'Data_out/Bubble_TwiHead_fulltable.xlsx', row.names = F)
#' Bubble_TwiHead_full$gene_short_name[which(Bubble_TwiHead_full$log2FoldChange < 0 & !is.na(Bubble_TwiHead_full$TF))]



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' comparision of Bubble and WThead; whether we see that twist is signigicantly different expressed;
res = results(object = DGE_bulkRNA, lfcThreshold = 1, contrast = c('phenotype', 'Bubble', 'WTHead'))

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

Bubble_WTHead = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
                       by.x = 'Gene', by.y = 'NV2Id', all.x = T)

#' refine based on padjust and logFC
Bubble_WTHead = Bubble_WTHead[which(abs(Bubble_WTHead$log2FoldChange) > 1 & Bubble_WTHead$padj < 0.01),, drop = F]
Bubble_WTHead = Bubble_WTHead[, -grep('^Twihead*', colnames(Bubble_WTHead))]
Bubble_WTHead = Bubble_WTHead[, -grep('Twi4d*', colnames(Bubble_WTHead))]
Bubble_WTHead = Bubble_WTHead[, -grep('WT4d*', colnames(Bubble_WTHead))]

#' add meta data to up-/down table
Bubble_WTHead_full = merge(x = Bubble_WTHead, y = NV2_annotation, 
                            by.x = 'Gene', by.y = 'NV2', all.x = T)
write.xlsx(x = Bubble_WTHead_full, file = 'Data_out/Bubble_WTHead_fulltable.xlsx', row.names = F)
#' Bubble_WTHead_full$gene_short_name[which(Bubble_WTHead_full$log2FoldChange < 0 & !is.na(Bubble_WTHead_full$TF))]




## Differential gene expression: 
#' DESeq function consits of the following steps:
#' estimation of size factors (controlling for differences in the sequencing depth of the samples);
#' estimation of dispersion values for each gene
#' fitting a generalized linear model (negative binomial model)
#' Making a different Metadata column (comparing all mutants to wildtype animals)
coldata = DataFrame(genotype = factor(c(rep('TwistMutant', 9), rep('WT', 6))))

# creating DESeq2 object:
bulkRNA_object_genotype = DESeqDataSetFromMatrix(countData = bulkRNA_modi, 
                                                  colData = coldata,
                                                  design = ~genotype)

DGE_bulkRNA_genotype = DESeq(bulkRNA_object_genotype)


#' extract the results
res = results(object = DGE_bulkRNA_genotype, lfcThreshold = 1, contrast = c('genotype', 'TwistMutant', 'WT'))

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA_genotype, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

TwiMutants_WT = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
                      by.x = 'Gene', by.y = 'NV2Id', all.x = T)

#' refine based on padjust and logFC
TwiMutants_WT = TwiMutants_WT[which(abs(TwiMutants_WT$log2FoldChange) > 1 & TwiMutants_WT$padj < 0.01),, drop = F]

#' add meta data to up-/down table
TwiMutants_WT_full = merge(x = TwiMutants_WT, y = NV2_annotation, 
                           by.x = 'Gene', by.y = 'NV2', all.x = T)
write.xlsx(x = TwiMutants_WT_full, file = 'Data_out/TwiMutants_WT_fulltable.xlsx', row.names = F)
TwiMutants_WT_full$gene_short_name[which(TwiMutants_WT_full$log2FoldChange < 0 & !is.na(TwiMutants_WT_full$TF))]




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Co-expression analysis with NvTwist
#' looking into MOX genes; NvMOXA-D
#' MOX genes are co-expressed with Twist
MOX_genes = NV2_annotation$NV2[grep(pattern = 'NvMOX.*', NV2_annotation$gene_short_name, ignore.case = T)]
transformed_data = log(counts(bulkRNA_object_phenotype, normalize = T))

MOX_frame = as.data.frame(transformed_data[row.names(transformed_data) %in% c('NV2.10864', MOX_genes), ])
MOX_frame$condition = c('Twist', rep('MOX_genes', 4))


#' TwiHead _ vs _ WTHead
mox_twihead_wthead = data.frame(condition = c(rep('TwiHead_WTHead', 30),
                                              rep('Bubble_TwiHead', 30),
                                              rep('Bubble_WTHead', 30)),
                                group = c(c(rep('Twist', 6), rep('MOX_genes', 24)),
                                               c(rep('Twist', 6), rep('MOX_genes', 24)),
                                          c(rep('Twist', 6), rep('MOX_genes', 24))),
                                value = c(7.139940,  6.964653,  6.947383, 5.999216, 5.994176, 5.911352,
                                          5.829818,  5.830702,  5.814066, 5.538977,  5.530957,  5.497548, 4.724877,  4.875909,  4.617815,
                                          4.762264,  4.280618,  4.377430, 6.039762, 5.977537, 6.023941, 5.893855, 5.900844, 5.866312,
                                          5.612539, 5.317396, 5.344587, 4.592015, 4.711260, 4.674014,
                                          6.432408, 6.363824, 6.367006, 7.139940,  6.964653,  6.947383,
                                          4.490474, 4.314435, 4.696607, 4.570517, 3.854903, 4.360135, 3.997997, 3.972686, 3.666988, 3.104179, 2.979434, 2.527553,
                                          5.829818,  5.830702,  5.814066, 5.538977,  5.530957,  5.497548, 4.724877,  4.875909,  4.617815, 4.762264,  4.280618,  4.377430,
                                          6.432408, 6.363824, 6.367006, 5.999216, 5.994176, 5.911352,
                                          4.490474, 4.314435, 4.696607, 4.570517, 3.854903, 4.360135, 3.997997, 3.972686, 3.666988, 3.104179, 2.979434, 2.527553,
                                          6.039762, 5.977537, 6.023941, 5.893855, 5.900844, 5.866312,
                                          5.612539, 5.317396, 5.344587, 4.592015, 4.711260, 4.674014))

mox_twihead_wthead$group = factor(mox_twihead_wthead$group, levels = c('Twist', 'MOX_genes'))

ggplot(mox_twihead_wthead, aes(x = group, y = value)) + 
  geom_jitter(width = 0.1) +
  facet_grid(~condition) +
  theme(aspect.ratio = 1) + theme_cowplot(font_size = 12) +
  theme(panel.border = element_rect(fill = NA, color = 'black', size = 1.5)) +
  labs(y = 'log-normalized counts', x = '')





#' Nematocytes marker genes;
Sox2 = NV2.18482
Ncol3 = NV2.10685
Carboxypeptidase (tbd).


Paraxis, tbx4/5? 
  
  #'NOTCH like members
  Serrate, Jagged?
  
  



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Misc
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



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Looking into contrasts; results(DGE_bulkRNA)
resultsNames(DGE_bulkRNA)

for(i in 2:length(resultsNames(DGE_bulkRNA))){
  print(resultsNames(DGE_bulkRNA)[i])
  var1 = strsplit(resultsNames(DGE_bulkRNA)[i], '_')[[1]][2]
  var2 = strsplit(resultsNames(DGE_bulkRNA)[i], '_')[[1]][4]
  
  #' running contrasts;
  res = results(DGE_bulkRNA,
                contrast = c('group', var1, var2),
                lfcThreshold = 1)
  res_out = merge(as.data.frame(res),
                  as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)),
                  by = 'row.names', sort = F)
  
  names(res_out)[1] = 'Gene'
  res_out = res_out[which(res_out$padj <= 0.01 & abs(res_out$log2FoldChange) > 1),, drop = F]
  write.table(x = res_out, file = paste0(var1, '_', var2, '_comparision.txt'), sep = '\t', quote = F, row.names = F)
}





