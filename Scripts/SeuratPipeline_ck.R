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
data.mutant = Read10X(data.dir = 'scRNA_tissue.pharynx.twistmutant/')
data.control = Read10X(data.dir = 'scRNA_tissue.pharynx.twistcontrol/')

lib1 = CreateSeuratObject(counts = data.mutant, project = 'tissue.pharynx.mutant')
lib2 = CreateSeuratObject(counts = data.control, project = 'tissue.pharynx.control')

#' calculate the fraction of mitochondrial involvment; extensive reads mapping to mitochondrial features
#' may indicate a dieing (apoptosis) cell;
mitochondria_IDs = NV2_genes$NV2[which(mitochondria %in% NV2_genes$gene_short_name)]
lib1[['percent.mt']] = PercentageFeatureSet(object = lib1, features = mitochondria_IDs)
lib2[['percent.mt']] = PercentageFeatureSet(object = lib2, features = mitochondria_IDs)

#' add library information
levels(lib1@meta.data$orig.ident) = 'Pharynx_Mutant'
levels(lib2@meta.data$orig.ident) = 'Pharynx_Control'

#' subset data based on some QC-metrics;
#' cells with < 200 features (genes) are filtered
#' cells with < 20,000 reads (nCountRNA) are filtered
#' cells where the percentage mitochondrial genes > 5.99 are filtered

# VlnPlot(object = lib1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
# VlnPlot(object = lib2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
lib1 = subset(x = lib1, subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 6)
lib2 = subset(x = lib2, subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 6)
  
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
  
