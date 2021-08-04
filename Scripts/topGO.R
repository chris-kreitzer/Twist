### load library
library(topGO)

### load gene -> GO terms
geneID2GO<-readMappings(file="full_genes.go.tsv")

### get a list of genenames
geneNames <- names(geneID2GO)

### read the list of genes in each of the files
mat<-read.table("mat.tsv", header=F)
myc <- read.table("myc.tsv", header=F)
sphase <- read.table("sphase.tsv", header=F)
mphase <- read.table("mphase.tsv", header=F)

### generate a binary list with the position of the selected genes in the gene universe
matList <- factor(as.integer(geneNames %in% mat$V1 ))
mycList <- factor(as.integer(geneNames %in% myc$V1 ))
sList <- factor(as.integer(geneNames %in% sphase$V1 ))
mList <- factor(as.integer(geneNames %in% mphase$V1 ))

### add the gene names
names(matList) <- geneNames
names(mycList) <- geneNames
names(sList) <- geneNames
names(mList) <- geneNames

###  MOLECULAR FUNCTION ENRICHMENT

###create the GOData object
matGO <- new("topGOdata", ontology = "MF", allGenes = matList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
mycGO <- new("topGOdata", ontology = "MF", allGenes = mycList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
sGO <- new("topGOdata", ontology = "MF", allGenes = sList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
mGO <- new("topGOdata", ontology = "MF", allGenes = mList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

### classic fisher test
matFis<- runTest(matGO, algorithm = "classic", statistic = "fisher")
matRes<-GenTable(matGO, classic=matFis,  ranksOf = "classic", topNodes = 20)
write.table(matRes, file="mat.MF.GO.sig.tsv", row.names=F, col.names=T)
mycFis<- runTest(mycGO, algorithm = "classic", statistic = "fisher")
mycRes<-GenTable(mycGO, classic=mycFis,  ranksOf = "classic", topNodes = 10)
write.table(mycRes, file="myc.MF.GO.sig.tsv", row.names=F, col.names=T)
mFis<- runTest(mGO, algorithm = "classic", statistic = "fisher")
mRes<-GenTable(mGO, classic=mFis,  ranksOf = "classic", topNodes = 20)
write.table(mRes, file="mphase.MF.GO.sig.tsv", row.names=F, col.names=T)
sFis<- runTest(sGO, algorithm = "classic", statistic = "fisher")
sRes<-GenTable(sGO, classic=sFis,  ranksOf = "classic", topNodes = 25)
write.table(sRes, file="sphase.MF.GO.sig.tsv", row.names=F, col.names=T)

### BIOLOGICAL PROCESS ENRICHMENT
matGO <- new("topGOdata", ontology = "BP", allGenes = matList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
matFis<- runTest(matGO, algorithm = "classic", statistic = "fisher")
matRes<-GenTable(matGO, classic=matFis,  ranksOf = "classic", topNodes = 100)
write.table(matRes, file="mat.BP.GO.sig.tsv", row.names=F, col.names=T)

mycGO <- new("topGOdata", ontology = "BP", allGenes = mycList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
mycFis<- runTest(mycGO, algorithm = "classic", statistic = "fisher")
mycRes<-GenTable(mycGO, classic=mycFis,  ranksOf = "classic", topNodes = 60)
write.table(mycRes, file="myc.BP.GO.sig.tsv", row.names=F, col.names=T)

mGO <- new("topGOdata", ontology = "BP", allGenes = mList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
mFis<- runTest(mGO, algorithm = "classic", statistic = "fisher")
mRes<-GenTable(mGO, classic=mFis,  ranksOf = "classic", topNodes = 250)
write.table(mRes, file="mphase.BP.GO.sig.tsv", row.names=F, col.names=T)

sGO <- new("topGOdata", ontology = "BP", allGenes = sList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
sFis<- runTest(sGO, algorithm = "classic", statistic = "fisher")
sRes<-GenTable(sGO, classic=sFis,  ranksOf = "classic", topNodes = 160)
write.table(sRes, file="sphase.BP.GO.sig.tsv", row.names=F, col.names=T)

### CELLULAR COMPONENT
matGO <- new("topGOdata", ontology = "CC", allGenes = matList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
matFis<- runTest(matGO, algorithm = "classic", statistic = "fisher")
matRes<-GenTable(matGO, classic=matFis,  ranksOf = "classic", topNodes = 10)
write.table(matRes, file="mat.CC.GO.sig.tsv", row.names=F, col.names=T)
mycGO <- new("topGOdata", ontology = "CC", allGenes = mycList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
mycFis<- runTest(mycGO, algorithm = "classic", statistic = "fisher")
mycRes<-GenTable(mycGO, classic=mycFis,  ranksOf = "classic", topNodes = 10)
write.table(mycRes, file="myc.CC.GO.sig.tsv", row.names=F, col.names=T)
mGO <- new("topGOdata", ontology = "CC", allGenes = mList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
mFis<- runTest(mGO, algorithm = "classic", statistic = "fisher")
mRes<-GenTable(mGO, classic=mFis,  ranksOf = "classic", topNodes = 40)
write.table(mRes, file="mphase.CC.GO.sig.tsv", row.names=F, col.names=T)
sGO <- new("topGOdata", ontology = "CC", allGenes = sList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
sFis<- runTest(sGO, algorithm = "classic", statistic = "fisher")
sRes<-GenTable(sGO, classic=sFis,  ranksOf = "classic", topNodes = 30)
write.table(sRes, file="sphase.CC.GO.sig.tsv", row.names=F, col.names=T)
