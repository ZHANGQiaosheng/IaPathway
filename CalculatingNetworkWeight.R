rm(list = ls())
gc()
rm()
setwd("E:/PhD/Manuscript/DATA/BRCA/NEW/scale")
library(readr)  
######## 1:37 is  normal samples; 38:70 is tumor samples
ExpMatrix <- read.delim("E:/PhD/Manuscript/DATA/BRCA/NEW/scale/brca.exp.txt", row.names=1, stringsAsFactors=FALSE)
MethMatrix <- read.csv("E:/PhD/Manuscript/DATA/BRCA/NEW/scale/BRCA.Meth.txt", row.names=1, sep="", stringsAsFactors=FALSE)
gene_Network <- read.delim("E:/PhD/Manuscript/DATA/BRCA/NEW/network/newnetwork.txt", stringsAsFactors=FALSE)
kcorr <- numeric(0)
scorr <- numeric(0)
for (i in 1:nrow(gene_Network)){
  logiclaTag <- !is.na(match(gene_Network$first[i], MethMatrix$Gene_Symbol)) && !is.na(match(gene_Network$first[i], row.names(ExpMatrix))) && !is.na(match(gene_Network$second[i], MethMatrix$Gene_Symbol)) && !is.na(match(gene_Network$second[i], row.names(ExpMatrix)))
  if(logiclaTag){
    geneMethNo <- which(as.character(gene_Network$first[i]) == MethMatrix$Gene_Symbol)
    geneMethMatrix <- MethMatrix[geneMethNo, 1:37]   ## number normal sample 
    geneMethMatrix <- t(geneMethMatrix)
    myPCA <- prcomp(geneMethMatrix, scale. = F, center = F) # Perform PCA
    vars <- apply(myPCA$x, 2, var)  
    props <- vars / sum(vars)
    Threshold <- min(which(cumsum(props)>0.95))
    geneMethMatrix95 <- myPCA$x[,1:Threshold]
    geneExpNo <- which(as.character(gene_Network$first[i]) == row.names(ExpMatrix)) # EXpression
    geneExpMatrix <- ExpMatrix[geneExpNo, 1:37]
    geneExpMatrix <- t(geneExpMatrix)
    mat1 <- cbind(geneMethMatrix95, geneExpMatrix)             # merge
    # mat1 <- scale(mat1)                                       # mat1
    geneMethNo <- which(as.character(gene_Network$second[i]) == MethMatrix$Gene_Symbol)
    geneMethMatrix <- MethMatrix[geneMethNo, 1:37]   ## number normal sample 
    geneMethMatrix <- t(geneMethMatrix)
    myPCA <- prcomp(geneMethMatrix, scale. = F, center = F) # Perform PCA
    vars <- apply(myPCA$x, 2, var)  
    props <- vars / sum(vars)
    Threshold <- min(which(cumsum(props)>0.95))
    geneMethMatrix95 <- myPCA$x[,1:Threshold]
    geneExpNo <- which(as.character(gene_Network$second[i]) == row.names(ExpMatrix)) # EXpression
    geneExpMatrix <- ExpMatrix[geneExpNo, 1:37]
    geneExpMatrix <- t(geneExpMatrix)
    mat2 <- cbind(geneMethMatrix95, geneExpMatrix)   
    # mat2 <- scale(mat2)                                       # mat2
    kca <- kcca(mat1, mat2, ncomps = 1)
    sca <- CCA(mat1, mat2, typex="standard", typez="standard", K=1)
    # ca <- cancor(mat1,mat2, xcenter = FALSE, ycenter = FALSE)
    # corcoef.test(r=ca$cor, n = nrow(mat1) , p = ncol(mat1), q = ncol(mat2))
    kcorr[i] <- kca@kcor
    scorr[i] <- sca$cors
  }else{
    kcorr[i] <- 0
    scorr[i] <- 0 
  }
}
gene.Network <- data.frame(gene_Network, kcorr, scorr)
getwd()
setwd("E:/PhD/Manuscript/DATA/BRCA/NEW/network/test")
write.table(gene.Network,"Nor.Weight.gene.Network.txt", col.names = T, row.names = T, sep = "\t", quote = F)

# plot(density(kcorr, bw=0.2))
# plot(density(scorr, bw=0.2))
# hist(kcorr)
# hist(scorr)
# 
# 
# colnames(MethMatrix)
# substr(colnames(ExpMatrix), c(1), c(16))
# substr(colnames(MethMatrix), c(1), c(16))
# 
# match(substr(colnames(ExpMatrix), c(1), c(16)), substr(colnames(MethMatrix), c(1), c(16)))
# library("PMA")
# CCA(mat1, mat2, typex="standard", typez="standard", K=1)
# library("kernlab")
# kcca(mat1, mat2, ncomps = 1)
