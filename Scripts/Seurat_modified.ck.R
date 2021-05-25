## scRNA analysis using Seurat;
## This script is an update of Alison's original script (provided via DropBox)
## 
## Start: 05/24
## 
## chris kreitzer


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
system('sshpass -p "chris2340" scp -r kreitzer@vlogin1.csb.univie.ac.at:/scratch/kreitzer/Nematostella/results/map/cellranger/polypcontrol/outs/filtered_feature_bc_matrix /Users/chriskreitzer/Documents/Github/Twist')


#here are the three parts of the script, can run all or some each time
setup = T #loads the gene lists - can save this as '.RData' and avoid running each time
PRE_ANALYSIS = T #this loads the raw data and generates Seurat files for each library
analysis = T #this analyzes the merged dataset; can also use this for each library separately
summarize.data.plots = T #this plots the summary plots from analysis
find.DEGs = T#This calculates DEG lists from clustering in 'analysis'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Annotation for SEURAT analysis: here, we use a pretty old annotation (previous to cupcake);
# update annotations for NV2 model (Juan)

if (setup) {
  #load and update gene names...
  genes = read_excel("D:/Dropbox/Nematostella/data/Import4R/genesfile/genes_mOs.xls") 
  #this is the current version
  genes = as.data.frame(genes)
  annotations = read_excel("D:/Dropbox/Nematostella/documents/Genelists/NVeannotations_NoDash.xls")
  known_genes = vector("double", nrow(annotations))
  for (i in 1:nrow(annotations)){
    known_genes[i] = match(annotations[i,1],genes[,1])
    if (as.logical(annotations[i,2]==".")){
      genes[known_genes[i],2] <- annotations[i,1]
      annotations[i,2] = annotations[i,1] #this takes care of indexing problem later for GO terms
      } else {
        genes[known_genes[i],2] <- annotations[i,2]
        genes[known_genes[i],3] <- annotations[i,3]
        genes[known_genes[i],4] <- annotations[i,5]
        genes[known_genes[i],5] <- annotations[i,7]
      }
    }
  
  
  # load TFs
  TF_NVE2 = read_excel("D:/Dropbox/Nematostella/documents/Genelists/tfs.xlsx")
  # #convert TF list to "names" 
  known_TFs = vector("double", nrow(TF_NVE2))
  
  for (i in 1:nrow(TF_NVE2)){
    known_TFs[i] = match(TF_NVE2$GeneId[i],genes[,1])
  }
  
  TF_list = genes[known_TFs, ]
  rm(TF_NVE2)
  
  # generate some gene lists for filtering:
  cellcyclearrest.genes = grep(pattern = "cell cycle arrest", annotations$gene_ontology_pfam)
  arrest_genes = vector("double", length(cellcyclearrest.genes))
  
  for(i in 1:length(cellcyclearrest.genes)){
    arrest_genes[i] = match(annotations[cellcyclearrest.genes[i], 1], genes[, 1])
  }
  cellcylearrest = genes$gene_short_name[arrest_genes]
  rm(cellcyclearrest.genes)
  
  mito.genes = grep(pattern = "mitochondrial", annotations$annotation_notes)
  mito_genes = vector("double", length(mito.genes)) 
  
  for(i in 1:length(mito.genes)){
    mito_genes[i] = match(annotations[mito.genes[i], 1], genes[, 1])
  }
  
  mitochondria = genes$gene_short_name[mito_genes]
  rm(mito.genes)
  
  Ribosomal.genes = grep(pattern = "ribosomal", annotations$annotation_notes)
  ribosomal_genes = vector("double", length(Ribosomal.genes)) 
  
  for(i in 1:length(Ribosomal.genes)){
    ribosomal_genes[i] = match(annotations[Ribosomal.genes[i], 1], genes[, 1])
  }
  
  ribosomal = genes[ribosomal_genes, 2]
  rm(Ribosomal.genes)
  
  save.image(file = '.RData')
  #save workspace here and then just load '.RData' each time you start and skip the begining
} else {
    load ('.RData')
  }



## loading data directly from LifeScienceCluster:
#' using system command to do so;

if (PRE_ANALYSIS)
{
  #direct to the matrix files of interest here.
  raw.data1 <- Read10X(data.dir = getwd())
  raw.data2 <- Read10X(data.dir = 'Z:/sequencing/Alison/Nematostella/TwistSequencing/ControlTissue_12000x2')
  
  # # set the gene names to the annotations
  rownames(raw.data1) <- genes$gene_short_name
  rownames(raw.data2) <- genes$gene_short_name
  
  #calculate mitochondrial fraction  
  percent.mito1 <- Matrix::colSums(raw.data1[mito_genes, ])/Matrix::colSums(raw.data1)
  percent.mito2 <- Matrix::colSums(raw.data2[mito_genes, ])/Matrix::colSums(raw.data2)
  
  #generate Seurat object    
  lib1  <- CreateSeuratObject(counts = raw.data1, project = "TWISTtissuesMutant")
  lib2  <- CreateSeuratObject(counts = raw.data2, project = "TWISTtissuesControl")
  
  #add mitochondria information
  lib1[["percent.mt"]] <- PercentageFeatureSet(object = lib1, features = mito_genes)
  lib2[["percent.mt"]] <- PercentageFeatureSet(object = lib2, features = mito_genes)
  
  #add library information
  levels(lib1@meta.data$orig.ident) <- 'TwistMutantTissue'
  levels(lib2@meta.data$orig.ident) <- 'ControlTissue'
  
  #filter the cells by genes detected 
  VlnPlot(object = lib1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
  lib1 <- subset(x = lib1, subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 6)
  
  VlnPlot(object = lib2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
  lib2 <- subset(x = lib2, subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 6)
  
  #add library info to names for later identification
  lib1 <- RenameCells(lib1, add.cell.id = "TwistMutantTissue")
  lib2 <- RenameCells(lib2, add.cell.id = "TwistControlTissue")
  
  #Normalize
  lib1 <- NormalizeData(lib1, scale.factor = 5000)
  lib2 <- NormalizeData(lib2, scale.factor = 5000)
  
  #merge two libraries  
  data1 = merge(lib1, lib2,merge.data = T)
  
  lib1<-FindVariableFeatures(lib1)
  lib2<-FindVariableFeatures(lib2)
  
  #generate CCA merge: this has changed
  # dataCCA<-RunCCA(lib1,lib2)
  
  #save your individual libraries so you don't have to do this again
  TwistTissueMutant = lib1
  TwistTissueControl = lib2
  save (TwistTissueMutant, file = 'TwistMutantTissue.Robj')
  save (TwistTissueControl, file = 'TwistControlTissue.Robj')
  save (data1, file = 'TWISTTissueMerged.raw.Robj')
  
  #clean up the workspace  
  rm (raw.data1, raw.data2, lib1, lib2)
}else{load (file = 'TWISTTissueMerged.raw.Robj')}

if (analysis)
{
  
  #look at quality of the dataset:  
  VlnPlot(object = data1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
  
  #calculate variable genes
  data1 <- FindVariableFeatures(data1)#, selection.method = 'mean.var.plot',mean.cutoff = c(0.1,Inf),dispersion.cutoff = c(0.1,Inf)
  VariableFeaturePlot(data1)
  
  #Scale the data
  data1 <- ScaleData(data1, features = data1@assays$RNA@var.features,
                     split.by = 'orig.ident') #this is important for batch effects
  #there are other methods for integrating libaries; can be explored.
  
  # run different reduction algorythms:
  #PCA
  data1 <- RunPCA(data1, pcs.compute = 50)
  ElbowPlot(object = data1, ndims = 50)
  d=21 #choose appropriate dimensions from the graph
  
  #UMAP
  data1 <- RunUMAP(data1, dims = 1:d,
                   reduction.name ='umap',reduction.key ='umap', #this is default; can change name to save different maps
                   #the following parameters control your output.
                   n.neighbors = 10L, #how many similar cells do you expect
                   spread =1, #how big is the axis space: 0.2 for lienages; 1 for separation
                   min.dist = 0.1, #how close to plot each cell higher=spread
                   local.connectivity = 100)#overall how connected is the graph
  
  #check that output is not based only on library or information content:
  library.plot=
    DimPlot(data1,group.by = 'orig.ident' ) + NoAxes()+labs(title = 'Library | ID')
  
  feature.plot = FeaturePlot(data1,'nFeature_RNA', 
                             cols = c('lightgrey','darkred'))+NoAxes()+labs(title = 'nFeatures')
  
  library.plot+feature.plot
  
  FeaturePlot(data1,
              c('NvTwist'),#
              cols =c('lightgrey',rev(brewer.pal(11 , "Spectral" ))),
              order = T,)&NoAxes()&NoLegend()
  
  #clustering: this is semi-independent of the UMAP, but UMAP is useful for visualizing the clustering
  # this can be discusssed.
  data1 <- FindNeighbors(object = data1,reduction ="pca",dims = 1:d,
                         nn.method = 'annoy',  
                         #these parameters are same as UMAP uses. This is NOT seurat default
                         annoy.metric = 'cosine', 
                         k.param = 20)
  clust.cp = unique (c(cols25(25),glasbey(32),alphabet2(26),alphabet(26)))
  #here you can calculate different clustering resolutions:  
  res = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2) 
  p=NULL
  for (i in 1:10)
  {
    data1 <- FindClusters(object = data1,resolution = res[i],random.seed = 0)
    p[[i]]=DimPlot(data1, label = T,label.size = 5, repel = T,
                   cols = clust.cp)+NoLegend()+NoAxes()+
      labs(title = c(res[i],'Clusters | ID'))
  } 
  
  #can use clustree to image clustering stability:
  library(clustree)
  clustree(data1, prefix = "RNA_snn_res.")
  
  #choose your desired resolution:
  p[[1]]+p[[2]]+p[[4]]+p[[5]]
  
  data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0)
  
  #calculate relationship between clusters; DISCUSS
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = TRUE, dims = c(1:d))
  library.plot =DimPlot(data1, label = F,
                        group.by = 'orig.ident')+NoAxes()+
    labs(title = 'Library | ID')+NoLegend()
  cluster.plot =DimPlot(data1, label = T,label.size = 5, pt.size = 1,
                        shape.by = 'orig.ident',
                        repel = T,
                        cols = clust.cp)+NoAxes()+
    labs(title = 'Clusters | ID')+NoLegend()
  library.plot+cluster.plot

  
  #can also generate a 3D plot:
  use.3d = F
  if (use.3d)
  {
    data1 <- RunUMAP(data1, n.neighbors = 10L,spread = 0.2,min.dist = 0.1,seed.use = 0, 
                     metric = 'cosine',local.connectivity = 200,reduction.name = 'umap3d',reduction.key = 'umap3d',#60
                     dims = c(1:d), n.components = 3)
    
    d13=DimPlot(data1, dims = c(1,3), reduction = 'umap3d', group.by = 'orig.ident')+NoLegend()+NoAxes()
    d23=DimPlot(data1, dims = c(2,3), reduction = 'umap3d', group.by = 'orig.ident',order = levels(data1$orig.ident))+NoLegend()+NoAxes()#,cols = rev(LibCP)
    d13+d23
    
    #if 3 dimensions:
    UMAP_1 <- data1@reductions$umap3d@cell.embeddings[,1]
    
    UMAP_2 <- data1@reductions$umap3d@cell.embeddings[,2]
    
    UMAP_3 <- data1@reductions$umap3d@cell.embeddings[,3]
    
    library(scatterplot3d)
    library(viridis)
    #colour library
    color.library = as.numeric(data1@meta.data$orig.ident)
    for (i in 1:length(LibCP))
    {  color.library[color.library==i]<-LibCP[i]}
    
    #colour clusters
    color = as.numeric(data1@active.ident)
    cp <- clust.cp[1:length(levels(data1@active.ident))]
    for (i in 1:length(cp))
    {  color[color==i]<-cp[i]}
    
    #colour gene
    GOI  ='NvTwist'#"memOrange2"
    color.gene = round(as.numeric(log2(1+data1@assays$RNA@counts[GOI,]))) #work around to show expression.
    # cp <-c('lightgrey',rev(brewer.pal(2 , "Spectral" )))
    cp <-  colorRampPalette(c('light grey', "dark red"))(max(color.gene)+1)
    for (i in 0:length(cp)) #works only for descrete colours; not sure how to code a gradient.
    {  color.gene[color.gene==i]<-cp[i+1]}
    
    library(rgl)
    plot.color = color.gene #here as cluster (color) or library (color.library) or gene (color.gene)
    par3d(windowRect = c(20, 30, 800, 800)) #set window size
    #plot an interactive 3d plot:
    plot3d(x = UMAP_1, y = UMAP_2, z = UMAP_3, setLab = T, xlab = NULL,ylab = NULL,zlab = NULL,
           box = F, lwd = 0, axes=F, tick.marks = F,grid = F,
           col = plot.color,type = "p", size = 3)
    
  }  
  #save your finalized object for future analyses with a proper name:
  Twist.Merge = data1
  save (Twist.Merge, file = 'TWISTTissueMerged.Robj')
}else{load (file = 'TWISTTissueMerged.Robj')
  data1=Twist.Merge}

if (summarize.data.plots)
{
  #view output
  # set a maximally unique colour palette for clustering:
  clust.cp = unique (c(cols25(25),glasbey(32),alphabet2(26),alphabet(26)))
  LibCP = c('grey','black')#any two contrasting colours
  ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
  colnames(ids.cluster.library) = c('ID','Library','CellCount')
  
  #generate plots:
  
  cluster.plot =DimPlot(data1, label = F,label.size = 5, pt.size = 1,
                        shape.by = 'orig.ident',
                        repel = T,
                        cols = clust.cp)+NoAxes()+
    labs(title = 'Clusters | ID')+NoLegend()
  
  dist.clust2=ggplot(ids.cluster.library, aes(fill=ID, y= (CellCount), x=as.integer(Library))) + 
    geom_area(position="fill", stat="identity",alpha=0.4 , size=.5, colour="white") +
    geom_bar(position="fill", stat="identity", width = 0.5)+
    scale_fill_manual(values = clust.cp)+
    ggtitle("Distribution of cell types in time and space")
  
  #visualize your plots; can also do these individually as desired
  # library.plot+dist.lib+
  dist.clust2+cluster.plot
}
  
if (find.DEGs)
{
  all.markers_TF <- FindAllMarkers(data1, features = TF_list$gene_short_name,
                                   only.pos = TRUE, min.cells.gene = 3,
                                   return.thresh = 0.001, max.cells.per.ident = 500)
  
  all.markers_variable <- FindAllMarkers(data1,logfc.threshold = 0.6,
                                         features = data1@assays$RNA@var.features,
                                         return.thresh = 0.001,
                                         only.pos = TRUE, max.cells.per.ident = 500)
  
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
