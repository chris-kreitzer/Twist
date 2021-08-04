## GeneSet Enrichment analysis:
## Start analysis data: 04/08/2021
## 
## chris kreitzer
## working with TopGo!


## Libraries and Dependencies
BiocManager::install('topGO')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('ALL')
BiocManager::install("hgu95av2.db")

.rs.restartR()
library(topGO)
library(ALL)

data(ALL)
data(geneList)

affyLib = paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)


sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)
