#################### TCGA2STAT expression data R3.3.2
setwd("E:/PhD/Manuscript/DATA/BRCA/NEW")
library("TCGA2STAT")  
Sys.setenv(TAR="C:/cygwin64/bin/tar", R_GZIPCMD="C:/cygwin64/bin/gzip")
brca.exp <- getTCGA(disease="BRCA", data.type="mRNA_Array")
save(brca.exp, file = "F:/Data/TCGA2STAT/brca_exp.rda")
rm(list = ls())
gc()
load("F:/Data/TCGA2STAT/BRCA/Original/brca_exp.rda")
brca.exp.tum.norm <- TumorNormalMatch(brca.exp$dat)
brca.exp.bytype <- SampleSplit(brca.exp$dat)
dimnames(brca.exp.bytype$normal)[[2]]
dimnames(brca.exp.tum.norm$primary.tumor)[[2]]
substr(dimnames(brca.exp.bytype$normal)[[2]],c(1),c(12))
NormalNO <- match(dimnames(brca.exp.tum.norm$primary.tumor)[[2]],substr(dimnames(brca.exp.bytype$normal)[[2]],c(1),c(12)))
TumorNO <- match(dimnames(brca.exp.tum.norm$primary.tumor)[[2]],substr(dimnames(brca.exp.bytype$primary.tumor)[[2]],c(1),c(12)))
barcode_N <- dimnames(brca.exp.bytype$normal)[[2]][NormalNO]
barcode_T <- dimnames(brca.exp.bytype$primary.tumor)[[2]][TumorNO]
write.table(barcode_N,"barcode_N.txt",col.names = F,row.names = F,quote = F)
write.table(barcode_T,"barcode_T.txt",col.names = F,row.names = F,quote = F)
########### TCGAbiolinks  methylation data R3.5.1
barcode_N <- read.table("E:/PhD/Manuscript/DATA/BRCA/NEW/barcode_N.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
barcode_T <- read.table("E:/PhD/Manuscript/DATA/BRCA/NEW/barcode_T.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
setwd("F:/Data")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
barcode1 <- barcode_N$V1
barcode <- substr(barcode1,c(1),c(16))
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  legacy = TRUE,
  barcode = barcode
)
BRCA.Meth.nor <- GDCprepare(query,summarizedExperiment = FALSE)

barcode1 <- barcode_T$V1
barcode <- substr(barcode1,c(1),c(16))
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  legacy = TRUE,
  barcode = barcode
)
BRCA.Meth.tum <- GDCprepare(query,summarizedExperiment = FALSE)
BRCA.Meth.nor <- BRCA.Meth.nor[, -(2:3)]
BRCA.Meth.tum <- BRCA.Meth.tum[, -(1:3)]
BRCA.Meth <- cbind(BRCA.Meth.nor, BRCA.Meth.tum) # merge nor and tum
BRCA.Meth <- na.omit(BRCA.Meth)
BRCA.Meth <- BRCA.Meth[!(BRCA.Meth$Gene_Symbol==""), ]
library(stringr)                                   ## get first Symbol
for (i in 1:length(BRCA.Meth$Gene_Symbol)) {
  BRCA.Meth$Gene_Symbol[i] <- str_split(BRCA.Meth$Gene_Symbol[i],";")[[1]][1]
}
Gene_Symbol <- BRCA.Meth$Gene_Symbol
getwd()
write.table(BRCA.Meth,"BRCA.Meth.txt", row.names = F, col.names = T, sep = "\t", quote = F)
BRCA.Meth$Gene_Symbol[1:10]
substr(colnames(BRCA.Meth.nor)[2:length(colnames(BRCA.Meth.nor))],c(1),c(16)) 
colnames(BRCA.Meth)

Gene_Symbol <- BRCA.Meth$Gene_Symbol
BRCA.Meth <- BRCA.Meth[, -(1)]
BRCA.Meth <- t(BRCA.Meth)
BRCA.Meth <-  scale(BRCA.Meth)
BRCA.Meth <- t(BRCA.Meth)
BRCA.Meth <- round(BRCA.Meth,4)
BRCA.Meth <- as.data.frame(BRCA.Meth)
BRCA.Meth$Gene_Symbol <- Gene_Symbol
write.table(BRCA.Meth,"./scale/BRCA.Meth.txt")
brca.exp <- brca.exp$dat
brca.exp <- as.data.frame(brca.exp)
colnames(brca.exp)
barcode1 <- barcode_N$V1
barcode <- substr(barcode1,c(1),c(16))
brca.exp.nor.NO <- match(barcode[match(substr(colnames(BRCA.Meth.nor)[2:length(colnames(BRCA.Meth.nor))],c(1),c(16)), barcode)
] ,substr(colnames(brca.exp),c(1),c(16)))
brca.exp.nor <- brca.exp[,brca.exp.nor.NO]
barcode1 <- barcode_T$V1
barcode <- substr(barcode1,c(1),c(16))
brca.exp.tum.NO <- match(barcode[match(substr(colnames(BRCA.Meth.tum)[1:length(colnames(BRCA.Meth.tum))],c(1),c(16)), barcode)
                                 ] ,substr(colnames(brca.exp),c(1),c(16)))
brca.exp.tum <- brca.exp[,brca.exp.tum.NO]

brca.exp <- cbind(brca.exp.nor, brca.exp.tum)
write.table(brca.exp,"brca.exp.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# rowMean replace NA  Expression
ind <- which(is.na(brca.exp.nor), arr.ind=TRUE)
brca.exp.nor[ind] <- rowMeans(brca.exp.nor, na.rm = TRUE)[ind[,1]]
ind <- which(is.na(brca.exp.tum), arr.ind=TRUE)
brca.exp.tum[ind] <- rowMeans(brca.exp.tum, na.rm = TRUE)[ind[,1]]
brca.exp <- cbind(brca.exp.nor, brca.exp.tum)
write.table(brca.exp,"./MeanreplaceNA/brca.exp.txt", row.names = T, col.names = T, sep = "\t", quote = F)
###
brca.exp <- t(brca.exp)
brca.exp <-  scale(brca.exp)
brca.exp <- t(brca.exp)
brca.exp <- as.data.frame(brca.exp)
write.table(brca.exp,"./scale/brca.exp.txt", row.names = T, col.names = T, sep = "\t", quote = F)


getwd()
rm(BRCA.Meth)


match(substr(colnames(BRCA.Meth.nor)[2:length(colnames(BRCA.Meth.nor))],c(1),c(16)), barcode)
match(substr(colnames(BRCA.Meth.tum)[1:length(colnames(BRCA.Meth.tum))],c(1),c(16)), barcode)
rm(list = ls())
gc()
