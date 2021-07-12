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

## Input and Processing
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
bulkRNA_object = DESeqDataSetFromMatrix(countData = bulkRNA_modi, 
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
vsd_bulkRNA = vst(bulkRNA_object, blind = FALSE)
rlog_bulkRNA = rlog(bulkRNA_object, blind = FALSE)

#' compare the normalization methods on twist;
bulkRNA_object = estimateSizeFactors(bulkRNA_object)

transformed_data = bind_rows(
  as_data_frame(log2(counts(bulkRNA_object, normalized = T)[, 1:15] + 1)) %>%
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
#' DESeq function consits of the following steps:
#' estimation of size factors (controlling for differences in the sequencing depth of the samples);
#' estimation of dispersion values for each gene
#' fitting a generalized linear model (negative binomial model)
DGE_bulkRNA = DESeq(bulkRNA_object)

#' get differential expression results;
#' Bubble vs WT4d
res = results(object = DGE_bulkRNA, lfcThreshold = 1, contrast = c('phenotype', 'Bubble', 'Twi4d'))
res = res[order(res$padj), ]

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

#' make annotations for the DE genes within the groups; focus on TF and muscle development;
ortho_table = read.csv(file = 'Data_out/NV2_orthologe.table.tsv', sep = '\t')
ortho_table = ortho_table[!duplicated(ortho_table$NV2Id) & !is.na(ortho_table$NV2Id), ]

#' one example: Twist 4d (mutant) early development vs Bubble (adult animal with mutation and phenotype)
Bubble_WT4d = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
          by.x = 'Gene', by.y = 'NV2Id', all.x = T)

#' refine based on padjust and logFC
Bubble_WT4d = Bubble_WT4d[which(abs(Bubble_WT4d$log2FoldChange) > 1 & Bubble_WT4d$padj < 0.01),, drop = F]
Bubble_WT4d = Bubble_WT4d[, -grep('^WThead*', colnames(Bubble_WT4d))]
Bubble_WT4d = Bubble_WT4d[, -grep('Twihead*', colnames(Bubble_WT4d))]
Bubble_WT4d = Bubble_WT4d[, -grep('Twi4d*', colnames(Bubble_WT4d))]

#' export results:
write.table(x = Bubble_WT4d, file = 'Data_out/Bubble_WT4d_dge_fulltable.txt', sep = '\t', row.names = F, quote = F)

#' just TF and muscle annotation
Bubble_WT4d_short = Bubble_WT4d[!is.na(Bubble_WT4d$TF) | !is.na(Bubble_WT4d$deM_TF), ]
write.table(x = Bubble_WT4d_short, file = 'Data_out/Bubble_WT4d_dge_filtered.txt', sep = '\t', row.names = F, quote = F)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' comparision of TwiHead and WT4d; whether we see that twist is signigicantly different expressed;
res = results(object = DGE_bulkRNA, lfcThreshold = 1, contrast = c('phenotype', 'TwiHead', 'WT4d'))
res = res[order(res$padj), ]

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

Twihead_WT4d = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
                    by.x = 'Gene', by.y = 'NV2Id', all.x = T)

#' refine based on padjust and logFC
Twihead_WT4d = Twihead_WT4d[which(abs(Twihead_WT4d$log2FoldChange) > 1 & Twihead_WT4d$padj < 0.01),, drop = F]
Twihead_WT4d = Twihead_WT4d[, -grep('^WThead*', colnames(Twihead_WT4d))]
Twihead_WT4d = Twihead_WT4d[, -grep('Twi4d*', colnames(Twihead_WT4d))]
Twihead_WT4d = Twihead_WT4d[, -grep('Bubble*', colnames(Twihead_WT4d))]

#' export results:
write.table(x = Twihead_WT4d, file = 'Data_out/Twihead_WT4d_dge_fulltable.txt', sep = '\t', row.names = F, quote = F)

#' just TF and muscle annotation
Twihead_WT4d_short = Twihead_WT4d[!is.na(Twihead_WT4d$TF) | !is.na(Twihead_WT4d$deM_TF), ]
write.table(x = Twihead_WT4d_short, file = 'Data_out/Twihead_WT4d_dge_filtered.txt', sep = '\t', row.names = F, quote = F)








#' Volcano plot with "significant" genes labeled
volcanoplot = function(res, lfcthresh = 1, sigthresh = 0.05, 
                       main = "Volcano Plot", 
                       legendpos = "bottomright", 
                       labelsig = TRUE, textcx = 1, ...) {
  
  with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = main, ...))
  with(subset(res, padj < sigthresh), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red", ...))
  with(subset(res, abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), pch = 20, col = "orange", ...))
  with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), pch = 20, col = "green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs = Gene, cex = textcx, ...))
  }
  legend(legendpos, xjust = 1, yjust = 1, 
         legend = c(paste("FDR<", sigthresh, sep = ""), 
                    paste("|LogFC|>", lfcthresh, sep = ""), "both"), 
         pch = 20, col = c("red", "orange", "green"))
}

png("diffexpr-volcanoplot.png", 1400, 1200, pointsize = 15)
volcanoplot(resdata, lfcthresh = 1, sigthresh = 0.05, textcx = .8, xlim = c(-2.3, 2))
dev.off()


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



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##' some other contrasts for DGE analysis:
res = results(object = DGE_bulkRNA, lfcThreshold = 1, contrast = c('phenotype', 'Twi4d', 'WT4d'))
summary(res)
res = res[order(res$padj), ]

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

#' make annotations for the DE genes within the groups; focus on TF and muscle development;
ortho_table = read.csv(file = 'Data_out/NV2_orthologe.table.tsv', sep = '\t')
ortho_table = ortho_table[!duplicated(ortho_table$NV2Id) & !is.na(ortho_table$NV2Id), ]

#' one example: Twist 4d (mutant) early development vs Bubble (adult animal with mutation and phenotype)
Twi4d_WT4d = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
                    by.x = 'Gene', by.y = 'NV2Id', all.x = T)
Twi4d_WT4d = Twi4d_WT4d[which(abs(Twi4d_WT4d$log2FoldChange) > 1 & Twi4d_WT4d$padj < 0.01),, drop = F]

#' refine based on padjust and logFC
Twi4d_WT4d = Twi4d_WT4d[, -grep('^WThead*', colnames(Twi4d_WT4d))]
Twi4d_WT4d = Twi4d_WT4d[, -grep('Twihead*', colnames(Twi4d_WT4d))]
Twi4d_WT4d = Twi4d_WT4d[, -grep('Bubble*', colnames(Twi4d_WT4d))]

write.xlsx2(Twi4d_WT4d, file = 'Data_out/Twi4d_WT4d.xlsx', sheetName = "DGE_Twi4d.WT4d",
            col.names = TRUE, row.names = F, append = FALSE)




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##' some other contrasts for DGE analysis:
res = results(object = DGE_bulkRNA, lfcThreshold = 1, contrast = c('phenotype', 'Bubble', 'TwiHead'))
summary(res)
res = res[order(res$padj), ]

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

#' make annotations for the DE genes within the groups; focus on TF and muscle development;
ortho_table = read.csv(file = 'Data_out/NV2_orthologe.table.tsv', sep = '\t')
ortho_table = ortho_table[!duplicated(ortho_table$NV2Id) & !is.na(ortho_table$NV2Id), ]

#' one example: Twist 4d (mutant) early development vs Bubble (adult animal with mutation and phenotype)
Bubble_Twihead = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
                   by.x = 'Gene', by.y = 'NV2Id', all.x = T)
Bubble_Twihead = Bubble_Twihead[which(abs(Bubble_Twihead$log2FoldChange) > 1 & Bubble_Twihead$padj < 0.01),, drop = F]

#' refine based on padjust and logFC
Bubble_Twihead = Bubble_Twihead[, -grep('^WThead*', colnames(Bubble_Twihead))]
Bubble_Twihead = Bubble_Twihead[, -grep('Twi4d*', colnames(Bubble_Twihead))]
Bubble_Twihead = Bubble_Twihead[, -grep('WT4d*', colnames(Bubble_Twihead))]

write.xlsx2(Bubble_Twihead, file = 'Data_out/Bubble_Twihead.xlsx', sheetName = "DGE_Bubble.Twihead",
            col.names = TRUE, row.names = F, append = FALSE)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##' some other contrasts for DGE analysis:
res = results(object = DGE_bulkRNA, lfcThreshold = 1, contrast = c('phenotype', 'Bubble', 'WTHead'))
summary(res)
res = res[order(res$padj), ]

#' Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(DGE_bulkRNA, normalized = TRUE)), 
                by = "row.names", sort = TRUE)
names(resdata)[1] = "Gene"

#' make annotations for the DE genes within the groups; focus on TF and muscle development;
ortho_table = read.csv(file = 'Data_out/NV2_orthologe.table.tsv', sep = '\t')
ortho_table = ortho_table[!duplicated(ortho_table$NV2Id) & !is.na(ortho_table$NV2Id), ]

#' one example: Twist 4d (mutant) early development vs Bubble (adult animal with mutation and phenotype)
Bubble_WThead = merge(resdata, ortho_table[,c('NV2Id', 'TF', 'deM_TF', 'BLAST.Hit', 'Trinotate.Descr', 'Emapper.Annotation')],
                       by.x = 'Gene', by.y = 'NV2Id', all.x = T)
Bubble_WThead = Bubble_WThead[which(abs(Bubble_WThead$log2FoldChange) > 1 & Bubble_WThead$padj < 0.01),, drop = F]

#' refine based on padjust and logFC
Bubble_WThead = Bubble_WThead[, -grep('Twihea*', colnames(Bubble_WThead))]
Bubble_WThead = Bubble_WThead[, -grep('Twi4d*', colnames(Bubble_WThead))]
Bubble_WThead = Bubble_WThead[, -grep('WT4d*', colnames(Bubble_WThead))]

write.xlsx2(Bubble_WThead, file = 'Data_out/Bubble_WThead.xlsx', sheetName = "DGE_Bubble.WThead",
            col.names = TRUE, row.names = F, append = FALSE)









