## The generell structure of the Seurat scRNA pipeline: 
## 
## Annotations: NV2 annotations (we need gene_short_names, 
## we need gene_ontology_pfam (cell-cycle, mito-genes, ribosomal genes)) - speak with Juan
## I made already the update (06/14/2021): load(NV2_annotations)
## 
## PRE-analyis:
##    exclude cells with too much mitochondria signal
##    normalization
##    several libraries
##    
## ANALYSIS Part:

rm(list = ls())
.rs.restartR()
setwd('~/Documents/GitHub/Twist/')



## Libraries
library(easypackages)
library(Seurat)
libraries("Seurat", 
          "Matrix", 
          "readxl",
          "RColorBrewer",
          'pals',
          'patchwork',
          'ggplot2')


## Input; Data load via system command:
## The system-command must be done just once (to load the data from the server to local machine)
load('NV2_annotation')
system('sshpass -p "chris2340" scp -r kreitzer@vlogin1.csb.univie.ac.at:/scratch/kreitzer/Nematostella/results/map/cellranger/polypcontrol/outs/filtered_feature_bc_matrix /Users/chriskreitzer/Documents/Github/Twist')
system('sshpass -p "chris2340" scp -r kreitzer@vlogin1.csb.univie.ac.at:/scratch/kreitzer/Nematostella/results/map/cellranger/tissue_pharynx_twistcontrol/outs/filtered_feature_bc_matrix /Users/chriskreitzer/Documents/Github/Twist')
system('sshpass -p "chris2340" scp -r kreitzer@vlogin1.csb.univie.ac.at:/scratch/kreitzer/Nematostella/results/map/cellranger/tissue_pharynx_twistmutant/outs/filtered_feature_bc_matrix /Users/chriskreitzer/Documents/Github/Twist')
system('sshpass -p "chris2340" scp -r kreitzer@vlogin1.csb.univie.ac.at:/scratch/kreitzer/Nematostella/results/map/cellranger/twistmutant_primarypolyp/outs/filtered_feature_bc_matrix /Users/chriskreitzer/Documents/Github/Twist')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Pre-analysis: 
#' Read in the data and create Seurat Object
data.mutant = Read10X(data.dir = 'scRNA_tissue.pharynx.twistcontrol/')
data.control = Read10X(data.dir = 'scRNA_tissue.pharynx.twistmutant/')

#' assign gene names rather than NV2.IDs
lib1 = data.mutant
lib2 = data.control

for(i in 1:length(row.names(lib1))){
  if(row.names(lib1)[i] %in% NV2_genes$NV2){
    row.names(lib1)[i] = NV2_genes$gene_short_name[which(NV2_genes$NV2 == row.names(lib1)[i])]
  } else {
    next
  }
}

#' there are non-unique entries in row.names(lib1):
all.row.names = row.names(lib1)
non.unique = row.names(lib1)[which(duplicated(row.names(lib1)))]

for(i in 1:length(non.unique)){
  duplicated.rownames = row.names(lib1)[which(row.names(lib1) == non.unique[i])]
  row.names(lib1)[which(row.names(lib1) == non.unique[i])] = paste0(duplicated.rownames, '-', letters[1:length(duplicated.rownames)])
}

#' second library:
for(i in 1:length(row.names(lib2))){
  if(row.names(lib2)[i] %in% NV2_genes$NV2){
    row.names(lib2)[i] = NV2_genes$gene_short_name[which(NV2_genes$NV2 == row.names(lib2)[i])]
  } else {
    next
  }
}

for(i in 1:length(non.unique)){
  duplicated.rownames = row.names(lib2)[which(row.names(lib2) == non.unique[i])]
  row.names(lib2)[which(row.names(lib2) == non.unique[i])] = paste0(duplicated.rownames, '-', letters[1:length(duplicated.rownames)])
}


#' create Seurat object
lib1 = CreateSeuratObject(counts = lib1, project = 'tissue.pharynx.mutant')
lib2 = CreateSeuratObject(counts = lib2, project = 'tissue.pharynx.control')

#' calculate the fraction of mitochondrial involvment; extensive reads mapping to mitochondrial features
#' may indicate a dieing (apoptosis) cell;
mitochondria = intersect(mitochondria, row.names(lib1))
lib1[['percent.mt']] = PercentageFeatureSet(object = lib1, features = mitochondria)
lib2[['percent.mt']] = PercentageFeatureSet(object = lib2, features = mitochondria)

#' add library information
levels(lib1@meta.data$orig.ident) = 'Pharynx_Mutant'
levels(lib2@meta.data$orig.ident) = 'Pharynx_Control'

#' subset data based on some QC-metrics;
#' cells with < 200 features (genes) are filtered
#' cells with < 20,000 reads (nCountRNA) are filtered
#' cells where the percentage mitochondrial genes > 5.99 are filtered

VlnPlot(object = lib1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
VlnPlot(object = lib2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
lib1 = subset(x = lib1, subset = nFeature_RNA > 250 & nCount_RNA < 20000 & percent.mt < 6)
lib2 = subset(x = lib2, subset = nFeature_RNA > 250 & nCount_RNA < 20000 & percent.mt < 6)
  
#' add library info to names for later identification
lib1 = RenameCells(lib1, add.cell.id = "TwistMutantTissue")
lib2 = RenameCells(lib2, add.cell.id = "TwistControlTissue")
  
#' Normalize; the expression of every feature in every cell;
#' We employ a global-scaling normalization method “LogNormalize” that normalizes 
#' the feature expression measurements for each cell by the total expression, 
#' multiplies this by a scale factor (5,000 by default), and log-transforms the result. 
#' Normalized values are stored in lib[["RNA"]]@data
lib1 = NormalizeData(lib1, scale.factor = 5000)
lib2 = NormalizeData(lib2, scale.factor = 5000)
  
#' merge two libraries  
data1 = merge(lib1, lib2, merge.data = T)
lib1 = FindVariableFeatures(lib1)
lib2 = FindVariableFeatures(lib2)

# TwistTissueMutant = lib1
# TwistTissueControl = lib2
# save (TwistTissueMutant, file = 'TwistMutantTissue.Robj')
# save (TwistTissueControl, file = 'TwistControlTissue.Robj')
# save (data1, file = 'TWISTTissueMerged.raw.Robj')
  
plot1 = VariableFeaturePlot(lib1)
top10 = head(VariableFeatures(lib1), 20)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

#' patchwork approach (package to arrange multiple plots)
plot1 / plot2


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ANALYSIS Part:
#'
#' Scaling the data (expression)
#' We apply a linear transformation (‘scaling’) 
#' that is a standard pre-processing step prior to dimensional reduction techniques like PCA. 
#' The ScaleData function shifts the expression of each gene, so that the mean expression across cells is 0
#' Scales the expression of each gene, so that the variance across cells is 1


#' calculate variable genes:
data1 = FindVariableFeatures(data1)

#' scale the data
data1 = ScaleData(data1, features = data1@assays$RNA@var.features,
                   split.by = 'orig.ident') # this is important for batch effects


#' PCA
data1 = RunPCA(data1, pcs.compute = 50)
ElbowPlot(object = data1, ndims = 50)
d = 18 # choose appropriate dimensions from the graph

#UMAP
data1 = RunUMAP(data1, 
                dims = 1:d,
                n.neighbors = 10L, # how many similar cells do you expect
                spread = 1, # how big is the axis space: 0.2 for lienages; 1 for separation
                min.dist = 0.1, # how close to plot each cell higher=spread
                local.connectivity = 100) # overall how connected is the graph

#check that output is not based only on library or information content:
library.plot = DimPlot(data1, 
                       group.by = 'orig.ident',
                       reduction = 'umap',
                       label = F) + 
  labs(title = 'Library | ID') +
  NoAxes()


feature.plot = FeaturePlot(data1,
                           'nFeature_RNA', 
                           cols = c('lightgrey','darkred')) +
  NoAxes() + 
  labs(title = 'nFeatures')

library.plot+feature.plot


FeaturePlot(data1,
            c('NvTwist'),
            shape.by = NULL,
            #label = T,
            cols = c('lightgrey', 
                     rev(brewer.pal(11 , "Spectral"))), order = T, split.by = 'orig.ident') & 
  NoAxes()
  


#' Clustering to cells:
data1 = FindNeighbors(object = data1,
                      reduction = "pca",
                      dims = 1:d,
                      nn.method = 'annoy',  
                      annoy.metric = 'cosine', 
                      k.param = 20)

clust.cp = unique(c(cols25(25), glasbey(32), alphabet2(26), alphabet(26)))

#here you can calculate different clustering resolutions:  
res = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2) 
p = NULL
for(i in 1:10){
  data1 = FindClusters(object = data1,
                       resolution = res[i],
                       random.seed = 0)
  p[[i]] = DimPlot(data1, 
                   label = T,
                   label.size = 5, 
                   repel = T,
                   cols = clust.cp) + 
    NoLegend() + 
    NoAxes() + 
    labs(title = c(res[i], 'Clusters | ID'))
}

p[[8]]


#can use clustree to image clustering stability:
library(clustree)
clustree::clustree(data1, prefix = "RNA_snn_res.")

# choose your desired resolution:
p[[1]] + p[[2]] + p[[4]] + p[[5]]

data1 = FindClusters(object = data1,
                     resolution = 0.8,
                     random.seed = 0)


#calculate relationship between clusters; DISCUSS
data1 = BuildClusterTree(object = data1, 
                         reorder = TRUE,
                         reorder.numeric = TRUE, 
                         dims = c(1:d))

library.plot = DimPlot(data1, 
                       label = F,
                       group.by = 'orig.ident') + 
  NoAxes() +
  labs(title = 'Library | ID')


cluster.plot = DimPlot(data1, 
                       label = T,
                       label.size = 5, 
                       pt.size = 1,
                       shape.by = 'orig.ident',
                       repel = T,
                       cols = clust.cp) + 
  NoAxes() + 
  labs(title = 'Clusters | ID') + 
  NoLegend()

library.plot + cluster.plot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clust.cp = unique(c(cols25(25), glasbey(32), alphabet2(26), alphabet(26)))
LibCP = c('grey', 'black')

ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
colnames(ids.cluster.library) = c('ID', 'Library', 'CellCount')

#generate plots:
cluster.plot = DimPlot(data1, 
                       label = F, 
                       label.size = 5, 
                       pt.size = 1,
                       shape.by = 'orig.ident',
                       repel = T,
                       cols = clust.cp) +
  NoAxes() + 
  labs(title = 'Clusters | ID')

dist.clust2 = ggplot(ids.cluster.library, aes(x=as.integer(Library), y = CellCount, fill = ID)) + 
  geom_area(position = "stack", stat = "identity", alpha = 0.4 , size = .5, colour = "white") +
  geom_bar(position= "stack", stat = "identity", width = 0.5)+
  scale_fill_manual(values = clust.cp)+
  ggtitle("Distribution of cell types in time and space")

#visualize your plots; can also do these individually as desired
# library.plot+dist.lib+
dist.clust2 + cluster.plot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all.markers_TF = FindAllMarkers(data1, 
                                features = NV2_TF$gene_short_name,
                                only.pos = TRUE, 
                                min.cells.gene = 3,
                                return.thresh = 0.001, 
                                max.cells.per.ident = 500)
  
all.markers_variable = FindAllMarkers(data1,
                                      logfc.threshold = 0.6,
                                      features = data1@assays$RNA@var.features,
                                      return.thresh = 0.001,
                                      only.pos = TRUE, 
                                      max.cells.per.ident = 500)
  
  all.markers <- FindAllMarkers(data1,logfc.threshold = 0.6,
                                return.thresh = 0.001,
                                only.pos = F, max.cells.per.ident = 500)
  
  #generate a collated list of unique DE genes
  
  list = NULL
  for (i in 1:length(levels(data1@active.ident)))
  {
    x=all.markers_variable[as.numeric(all.markers_variable$cluster)==i,][1:min(10,length(which(as.numeric(all.markers_variable$cluster)==i))),7]
    if (is.na (x) ==F)
      list=c(list,x)
  }
  
  #Image the list
  DEG.variable= DotPlot(data1, features = unique(c(list)), 
                        scale.by='size' , col.min = 0, col.max = 3,  
                        cols = c('lightgrey','darkred')) + RotatedAxis() +
    FontSize(6,6)+NoLegend()+coord_flip()
  
  #also for the TF list
  list_TF = NULL
  for (i in 1:length(levels(data1@active.ident)))
  {
    x=all.markers_TF[as.numeric(all.markers_TF$cluster)==i,][1:min(10,length(which(as.numeric(all.markers_TF$cluster)==i))),7]
    if (is.na (x) ==F)
      list_TF=c(list_TF,x)
  }
  DEG.TFs=DotPlot(data1, features = unique(c(list_TF,'NvTwist')), 
                  scale.by='size' , col.min = 0, col.max = 3, 
                  cols = c('lightgrey','darkred')) + RotatedAxis() +
    FontSize(6,6)+NoLegend()+coord_flip()
  
  DEG.variable+DEG.TFs
  
}



