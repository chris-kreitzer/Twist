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

write.table(x = cluster20, file = 'Data_out/ConservedMarkers_CellCluster20_mutant.txt', row.names = T, sep = '\t')


#' finding markers which distinguish the both libraries.
markers = FindMarkers(data1,
                      min.pct = 0.45,
                      logfc.threshold = 1,
                      ident.1 = 'Pharynx_Mutant', 
                      ident.2 = 'Pharynx_Control',
                      group.by = 'orig.ident', 
                      subset.ident = '20')


#' features which are depleted and enriched in either mutant or wt animals;
#' when concentrating on cell-cluster 20
GOI.20 = row.names(markers)

plot_list = list()
for(i in unique(GOI.20)){
  plot = VlnPlot(data1, 
                 features = i,
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
  
  plot_list[[i]] = plot
  
}
  
upregulated = plot_list$NvNcol6 | plot_list$NVE8163 | plot_list$`RS19-like-a`
downregulated = plot_list$`FRIS-like5;FRIS-like7` | plot_list$`BAG4-like1` | plot_list$`ASGL1-like3-b`


upregulated / downregulated


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




## further:
#' look cell cluster identity; assign cell types
#' look into different features, across different conditions (mutant and wild-type)
#' work still on cluster 20 and 13
#' look if tentacle identity



## Identify conserved cell type markers
conserved.markers = FindConservedMarkers(data1, 
                                         ident.1 = 20, 
                                         grouping.var = "orig.ident", verbose = FALSE)
View(conserved.markers)

FeaturePlot(data1, features = c('NVE19039', 'FRIS-like4;FRIS-like6', 
                                'NVE10241-a', 'Nve-MELC4', 
                                'CSL3-like6', 'HE-like2'), 
            min.cutoff = "q9")




#' test
NVE8917 = VlnPlot(data1, features = c('TBB-like4;TBB4-like2'),
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


NVE8917















