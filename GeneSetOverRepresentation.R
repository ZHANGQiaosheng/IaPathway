rm(list = ls())
gc()
source("https://bioconductor.org/biocLite.R")
biocLite("HTSanalyzeR")
library("HTSanalyzeR")
library("limma")
eSet <- read.delim("E:/PhD/Manuscript/DATA/BRCA/NEW/scale/brca.exp.txt", row.names=1, stringsAsFactors=FALSE)
normals <- c(rep("Nor", 37), rep("Tur", 33))
design <- model.matrix(~0 + factor(normals))
colnames(design) <- c("Nor","Tur")
contrast.matrix <- makeContrasts(Nor-Tur, levels=design)
##step1
fit <- lmFit(eSet, design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
#tempOutput = topTable(fit2, coef=1,  n=Inf)
tempOutput = topTable(fit2, coef=1, p.value=0.01, n=Inf)
PathwayList <- list()
for (i in 1:281) {
  t1 <-
    paste(
      "E:/PhD/Manuscript/DATA/BRCA/Result/OldPathway/Original/pathwayinformation/pathwayNo/",
      i,
      sep = ''
    )
  path <-
    read.table(
      paste(t1, ".txt", sep = ''),
      quote = "\"",
      comment.char = "",
      as.is = TRUE
    )
  PathwayList[[i]] <- path$V1
}
pathwayID <- read.delim("E:/PhD/Manuscript/DATA/BRCA/Result/OldPathway/pathwayinformation/pathwayID.txt", stringsAsFactors=FALSE)
# pathwaynames <- list()
# for (i in 1:281){
#   pathwaynames[[i]]<-pathwayID$Symbel[i]
# }
names(PathwayList) <- pathwayID$Symbel
Result <- multiHyperGeoTest(PathwayList, row.names(eSet), row.names(tempOutput), minGeneSetSize = 1, pAdjustMethod = "BH", verbose = TRUE)
write.table(Result, "Result.txt", row.names = T, col.names = T, sep = "\t", quote = F)
