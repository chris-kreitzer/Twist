## Filtering INDEL variants called by bcftools;
## 
## 05/11/2021:



rm(list = ls())
.rs.restartR()
setwd('~/Documents/GitHub/Twist/')


## Libraries
library(stringr)


## Input
chr2 = read.csv('Data_out/chr2_filtered_intersected.txt', sep = '\t', header = F)
colnames(chr2)[1:5] = c('Chromsome', 'Position', 'ID', 'ref', 'alt')


## Processing:
FilterVariants = function(data.in){
  data.in = as.data.frame(data.in)
  colnames(data.in)[1:5] = c('Chromsome', 'Position', 'ID', 'ref', 'alt')
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
  
  # filter variants where more than 2 samples are heterozygous or non-genotyped;
  for(i in c(10:24)){
    data.single[, i] = substr(x = data.single[, i], start = 1, stop = 3)
  }
  
  data.single$het = NA
  data.single$het = apply(data.single[10:24], 1, function(x) length(which(x == '0/1')))
  
  data.single$keep = ifelse(is.na(data.single$filter), 'discard',
                            ifelse(data.single$filter == 'keep' & data.single$het > 2, 'discard', 'keep'))
  
  
  ## double, or more alt alleles:
  data.double = data.in[which(data.in$length > 1), ]
  data.double$INS = NA
  data.double$DEL = NA
  
  for(i in 1:nrow(data.double)){
    data.double$INS[i] = any(grepl(pattern = data.double$ref[i], x = unlist(strsplit(data.double$alt[i], split = ',', fixed = T))))
    data.double$DEL[i] = any(grepl(pattern = paste(unlist(strsplit(data.double$alt[i], split = ',', fixed = T)), collapse = '|'),
                                   x = data.double$ref[i]))
  }
  
  
  
  data.double$delta = max(abs(nchar(example$ref) - nchar(unlist(strsplit(example$alt, split = ',')))))
  data.double$keep = ifelse(data.double$match & data.double$delta > 3, 'keep', 'discard')
  
  # merge the two data frames
  data.out = rbind(data.single, data.double)
  
  return(data.out)
  
}


## Filtering heterozygous variants 
