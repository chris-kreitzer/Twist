## Additional downstream analysis for scRNA-seq in Pharynx-mutants and control animals

#' keep this plot for further analysis and visualization
plotB = DimPlot(data1, reduction = "umap", 
        group.by = "orig.ident", 
        cols = c('#CD9B06', '#4A265A'), 
        shuffle = T,
        label = F) +
  NoAxes()

plotA = DimPlot(data1, reduction = "umap") +
  NoAxes()

plotA + plotB

#' find differentially expressed features (cluster biomarkers)
#' Cell-Cluster 20, which is obviously different from the mutant and Control
cluster20 = FindMarkers(data1, 
                        ident.1 = 20, 
                        min.pct = 0.25, 
                        logfc.threshold = 1) #' conservative approach as logFC is set to 1


#' enriched features in cluster 20
#' here concentrating on NvNcol6
NvNcol6 = VlnPlot(data1, features = c('NvNcol6'),
                  split.by = 'orig.ident', 
                  split.plot = T) +
  RestoreLegend(position = "right") +
  labs(x = 'Cell-cluster identity') +
  scale_y_continuous(expand = c(0,0)) +
  theme(aspect.ratio = 0.5,
        axis.line = element_line(size = 0.9),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        legend.position = c(0.05, 0.9)) +
  scale_fill_manual(values = c('Pharynx_Control' = '#CD9B06', 
                                'Pharynx_Mutant' = '#4A265A')) +
  panel_border(color = 'black', size = 1.2)
  

NvNcol6

#' further upregulated gene
'NvNEP3-b'
#' downregulated gene
'ASGL1-like3-b', 'BAG4-like1'

# VlnPlot(data1, features = c('ASGL1-like3-b'),
#         split.by = 'orig.ident', 
#         split.plot = T) +
#   RestoreLegend(position = "right") +
#   labs(x = 'Cell-cluster identity') +
#   scale_y_continuous(expand = c(0,0)) +
#   theme(aspect.ratio = 0.5,
#         axis.line = element_line(size = 0.9),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16, face = 'bold'),
#         legend.position = c(0.05, 0.9)) +
#   scale_fill_manual(values = c('Pharynx_Control' = '#CD9B06', 
#                                'Pharynx_Mutant' = '#4A265A')) +
#   panel_border(color = 'black', size = 1.2)










