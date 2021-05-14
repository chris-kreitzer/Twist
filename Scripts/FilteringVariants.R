## Filtering INDEL variants called by bcftools;
## Pre-defined filtering steps, where we believe that variants will not effect protein function
## Concentrate on unique and multipe alternative genoypes respectively
## Filtering rules, as commented below;
##
## Start: 05/11/2021:
## Chris Kreitzer



rm(list = ls())
.rs.restartR()
setwd('~/Documents/GitHub/Twist/')


## Libraries
library(stringr)
library(tidyverse)

## Input
chr2 = read.csv('Data_out/chr2_filtered_intersected.txt', sep = '\t', header = F)



## Processing
FilterVariants = function(data.in){
  columns.keep = c('Chromosome', 'Position', 'ref', 'alt', 'INS', 'DEL', 'delta', 'het', 'keep', 'Annotation', 'gene_id')
  data.in = as.data.frame(data.in)
  colnames(data.in)[1:5] = c('Chromosome', 'Position', 'ID', 'ref', 'alt')
  annotation = grep('gene_id*', data.in)
  colnames(data.in)[annotation] = 'Annotation'
  data.in$length = lengths(strsplit(data.in$alt, split = ',', fixed = T))
  
  
  ## single alt alleles:
  data.single = data.in[which(data.in$length == 1), ]
  data.single$INS = NA
  data.single$DEL = NA
  
  for(i in 1:nrow(data.single)){
    data.single$INS[i] = grepl(pattern = data.single$ref[i], x = data.single$alt[i], fixed = T)
    data.single$DEL[i] = grepl(pattern = data.single$alt[i], x = data.single$ref[i], fixed = T)
  }
  
  # annotate variants which needs to be filtered
  data.single$delta = abs(nchar(data.single$ref) - nchar(data.single$alt))
  
  data.single$filter = NA
  ii.mismatch = which(data.single$INS == data.single$DEL)
  ii.match = which(data.single$INS != data.single$DEL & data.single$delta > 3)
  
  data.single[c(ii.mismatch, ii.match), 'filter'] = 'keep'
  
  data.single$repeatsREF = NA
  data.single$repeatsALT = NA
  for(i in 1:nrow(data.single)){
    data.single$repeatsREF[i] = sum(!!str_count(data.single$ref[i], LETTERS))
    data.single$repeatsALT[i] = sum(!!str_count(data.single$alt[i], LETTERS))
  }
  
  
  # filter variants where more than 2 samples are heterozygous or non-genotyped;
  start.samples = grep('GT:PL', data.single) + 1
  end.samples = start.samples + 15
  
  for(i in seq(start.samples, end.samples, 1)){
    data.single[, i] = substr(x = data.single[, i], start = 1, stop = 3)
  }
  
  data.single$het = NA
  data.single$het = apply(data.single[seq(start.samples, end.samples, 1)], 1, function(x) length(which(x == '0/1')))
  
  data.single$keep = ifelse(is.na(data.single$filter), 'discard',
                            ifelse(data.single$filter == 'keep' & data.single$het > 2, 'discard', 
                                   ifelse(data.single$repeatsREF < 3, 'discard', 'keep')))
  
  data.single$gene_id = as.character(lapply(strsplit(as.character(data.single$Annotation), split = ';'), 
                                            head, n = 1))
  
  data.single = data.single[!duplicated(data.single[,c('Position', 'gene_id')]), ]
  
  data.single = data.single[, colnames(data.single) %in% columns.keep]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # double, or more alternative genotypes
  data.double = data.in[which(data.in$length > 1), ]
  max_alternative_genotypes = max(data.in$length)
  
  data.multipe.genotypes = data.double %>%
    separate(alt, c(paste0('Genotype', seq(1, max_alternative_genotypes, 1))), 
             sep = ',',
             remove = F)
  
  # determine whether Insertion or Deletion: is ref-pattern found within alternative alleles
  # working on Insertions (ref shorter than alternative genotypes)
  data.multipe.genotypes$INS = NA
  
  for(i in 1:nrow(data.multipe.genotypes)){
    data.multipe.genotypes$INS[i] = any(grepl(pattern = data.multipe.genotypes$ref[i],
                                              x = data.multipe.genotypes[i, c(grep(pattern = 'Genotype*', x = colnames(data.multipe.genotypes)))]))
  }
  
  # working on deletions: alterative genotypes shorter than reference
  data.multipe.genotypes$DEL = NA
  genotypes = gsub(pattern = ',', replacement = '|', x = data.multipe.genotypes$alt)
  data.multipe.genotypes$alt = genotypes
  
  for(i in 1:nrow(data.multipe.genotypes)){
    data.multipe.genotypes$DEL[i] = grepl(data.multipe.genotypes$alt[i], data.multipe.genotypes$ref[i])
  }
  
  # determine INDEL extend
  data.multipe.genotypes$delta = NA
  
  for(i in 1:nrow(data.multipe.genotypes)){
    data.multipe.genotypes$delta[i] = max(abs(nchar(data.multipe.genotypes$ref[i]) - 
                                                      nchar(data.multipe.genotypes[i, c(grep(pattern = 'Genotype*', x = colnames(data.multipe.genotypes)))])), na.rm = T)
  }
  
  data.multipe.genotypes$repeatsREF = NA
  for(i in 1:nrow(data.multipe.genotypes)){
    data.multipe.genotypes$repeatsREF[i] = sum(!!str_count(data.multipe.genotypes$ref[i], LETTERS))
  }
  
  
  ## filter variants where more than 2 samples are heterozygous or non-genotyped
  start.samples = grep('GT:PL', data.multipe.genotypes) + 1
  end.samples = start.samples + 15
  
  for(i in seq(start.samples, end.samples, 1)){
    data.multipe.genotypes[, i] = substr(x = data.multipe.genotypes[, i], start = 1, stop = 3)
  }
  
  data.multipe.genotypes$het = NA
  data.multipe.genotypes$het = apply(data.multipe.genotypes[seq(start.samples, end.samples, 1)], 1, function(x) length(which(x == '0/1')))
  
  data.multipe.genotypes$keep = ifelse(data.multipe.genotypes$delta > 3 & 
                                         data.multipe.genotypes$het < 3 & 
                                         data.multipe.genotypes$repeatsREF > 2, 'keep', 'discard')
  
  data.multipe.genotypes$gene_id = as.character(lapply(strsplit(as.character(data.multipe.genotypes$Annotation), split = ';'), 
                                            head, n = 1))
  
  data.multipe.genotypes = data.multipe.genotypes[!duplicated(data.multipe.genotypes[,c('Position', 'gene_id')]), ]
  
  data.multipe.genotypes = data.multipe.genotypes[, colnames(data.multipe.genotypes) %in% columns.keep]
  
  
  # combined output
  Variants_out = rbind(data.single, data.multipe.genotypes)
  Variants_keep = Variants_out[which(Variants_out$keep == 'keep'), ]
  
  # attach the respective genotypes
  Variants_annotated = merge(Variants_keep, data.in,
                             by.x = c('Position'), by.y = c('Position'), all.x = T)
  Variants_annotated = Variants_annotated[!duplicated(Variants_annotated[,c('Position', 'gene_id')]), ]
  
  return(list(Variants = Variants_keep,
              Variants_GT = Variants_annotated))
  
}



x = FilterVariants(data.in = chr2)
write.table(x$Variants, file = 'Data_out/Filtered_Annotated_Variants_chr2.xlsx', sep = '\t', row.names = F, quote = F)
write.table(x$Variants_GT, file = 'Data_out/Filtered_Annotated_Variants_chr2_Genotype.xlsx', sep = '\t', row.names = F, quote = F)
