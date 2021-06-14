## Making new annotation table for Seurat's scRNA analysis
## Transferring NVE --> NV2 IDs
## 
## Start: 06/10/2021:
## 
## chris kreitzer:


#' Approach:
#' starting with Juan's confident map


setwd('~/Documents/GitHub/Twist/')
.rs.restartR()

## Input
load('Scripts/RData')
Vienna_confident = read.csv('Data_out/vienna_nv2.confident.map', sep = '\t', header = F)
colnames(Vienna_confident) = c('NVE', 'NV2')
orthotable_juan = read.csv('Data_out/NV2_orthologe.table.tsv', sep = '\t')



## Processing:
#' subset the orignal annotation table, to those NVE IDs with confident map
#' Starting with annotation table
out = data.frame()
for(i in 1:nrow(annotations)){
  confident = Vienna_confident[which(Vienna_confident$NVE == annotations$NVE[i]), ]
  
  if(nrow(confident) == 0){
    next
  } else if(nrow(confident) == 1){
    modi = annotations[i, ]
    modi$NV2 = confident$NV2
  } else if(nrow(confident) > 1 & nrow(confident) < 5){
    modi = annotations[i, ]
    modi = rbind(modi, modi[rep(1, nrow(confident) - 1), ])
    modi$NV2 = confident$NV2
  } else {
    next
  }
  rm(confident)
  out = rbind(out, modi)
}
  
updated_annotations = out
updated_annotations$NVE = updated_annotations$NV2
updated_annotations$NV2 = NULL
colnames(updated_annotations)[1] = 'NV2'

#' merge cells;
NV2_annotation = data.frame()
for(i in unique(updated_annotations$NV2)){
  data.subset = updated_annotations[which(updated_annotations$NV2 == i), ]
  if(nrow(data.subset) > 1){
    updated = data.subset[1, ]
    updated$gene_ontology_pfam = paste(data.subset$gene_ontology_pfam, collapse = ';')
    updated$gene_short_name = paste(data.subset$gene_short_name, collapse = ';')
  } else {
    updated = data.subset
  }
  NV2_annotation = rbind(NV2_annotation, updated)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' gene table
out = data.frame()
for(i in 1:nrow(genes)){
  confident = Vienna_confident[which(Vienna_confident$NVE == genes$NVE[i]), ]
  
  if(nrow(confident) == 0){
    next
  } else if(nrow(confident) == 1){
    modi = genes[i, ]
    modi$NV2 = confident$NV2
  } else if(nrow(confident) > 1 & nrow(confident) < 5){
    modi = genes[i, ]
    modi = rbind(modi, modi[rep(1, nrow(confident) - 1), ])
    modi$NV2 = confident$NV2
  } else {
    next
  }
  rm(confident)
  out = rbind(out, modi)
}

updated_genes = out
updated_genes$NVE = updated_genes$NV2
updated_genes$NV2 = NULL
colnames(updated_genes)[1] = 'NV2'

#' merge cells;
NV2_genes = data.frame()
for(i in unique(updated_genes$NV2)){
  data.subset = updated_genes[which(updated_genes$NV2 == i), ]
  if(nrow(data.subset) > 1){
    updated = data.subset[1, ]
    updated$gene_short_name = paste(data.subset$gene_short_name, collapse = ';')
    updated$annotation_notes = paste(data.subset$annotation_notes, collapse = ';')
  } else {
    updated = data.subset
  }
  NV2_genes = rbind(NV2_genes, updated)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' TF list
out = data.frame()
for(i in 1:nrow(TF_list)){
  confident = Vienna_confident[which(Vienna_confident$NVE == TF_list$NVE[i]), ]
  
  if(nrow(confident) == 0){
    next
  } else if(nrow(confident) == 1){
    modi = TF_list[i, ]
    modi$NV2 = confident$NV2
  } else if(nrow(confident) > 1 & nrow(confident) < 5){
    modi = TF_list[i, ]
    modi = rbind(modi, modi[rep(1, nrow(confident) - 1), ])
    modi$NV2 = confident$NV2
  } else {
    next
  }
  rm(confident)
  out = rbind(out, modi)
}

updated_TF = out
updated_TF$NVE = updated_TF$NV2
updated_TF$NV2 = NULL
colnames(updated_TF)[1] = 'NV2'

#' merge cells;
NV2_TF = data.frame()
for(i in unique(updated_TF$NV2)){
  data.subset = updated_TF[which(updated_TF$NV2 == i), ]
  if(nrow(data.subset) > 1){
    updated = data.subset[1, ]
    updated$gene_short_name = paste(data.subset$gene_short_name, collapse = ';')
    updated$annotation_notes = paste(data.subset$annotation_notes, collapse = ';')
  } else {
    updated = data.subset
  }
  NV2_TF = rbind(NV2_TF, updated)
}
colnames(NV2_TF)[3] = 'annotationOrthMCL'

#' complement with ID's from juan's orthotable
TF_Juan = orthotable_juan[!is.na(orthotable_juan$TF) | !is.na(orthotable_juan$deM_TF),, drop = F]

missing_transcriptionfactor = data.frame()
for(i in unique(setdiff(TF_Juan$NV2Id, NV2_TF$NV2))){
  if(is.na(i)){
    next
  } else {
    missing_TF = TF_Juan[which(TF_Juan$NV2Id == i), ]
    missing_TF = missing_TF[, c('NV2Id', 'Trinotate.Descr', 'Emapper.Name')]
    missing_transcriptionfactor = rbind(missing_transcriptionfactor, missing_TF)
    }
  }

colnames(missing_transcriptionfactor) = c('NV2', 'annotation_notes', 'gene_short_name')
missing_transcriptionfactor$annotationOrthMCL = NA
missing_transcriptionfactor$annotation_mod = NA

TF_updated = rbind(NV2_TF, missing_transcriptionfactor)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' specific gene lists
#' cell cycle arrest
cellcyclearrest.genes = grep(pattern = "cell cycle arrest", NV2_annotation$gene_ontology_pfam)
arrest_genes = vector("double", length(cellcyclearrest.genes))
for (i in 1:length(cellcyclearrest.genes)){
  arrest_genes[i] = match(NV2_annotation[cellcyclearrest.genes[i], 1], NV2_genes[, 1])
  }
cellcyclearrest = NV2_genes$gene_short_name[arrest_genes]

#' mitochondrial genes
mito.genes = grep(pattern = "mitochondrial", NV2_annotation$annotation_notes)
mito_genes = vector("double", length(mito.genes)) 

for (i in 1:length(mito.genes)){
  mito_genes[i] = match(NV2_annotation[mito.genes[i], 1], NV2_genes[, 1])
  }
mitochondria = NV2_genes$gene_short_name[mito_genes]

#' ribosomal genes:
Ribosomal.genes = grep(pattern = "ribosomal", NV2_annotation$annotation_notes)
ribosomal_genes = vector("double", length(Ribosomal.genes)) 

for (i in 1:length(Ribosomal.genes)){
  ribosomal_genes[i] = match(NV2_annotation[Ribosomal.genes[i], 1], NV2_genes[, 1])
  }
ribosomal = NV2_genes[ribosomal_genes, 2]


save(c('NV2_annotation', 'NV2_genes', 'NV2_TF', 'cellcyclearrest', 'mitochondria', 'ribosomal'),
           file = 'NV2_annotation')

save(NV2_annotation, NV2_genes, NV2_TF, cellcyclearrest, mitochondria, ribosomal,
           file = 'NV2_annotation')






































