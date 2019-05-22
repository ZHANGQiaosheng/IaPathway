setwd("E:/PhD/Manuscript/DATA/BRCA/NEW/Result")
ExPathwayList <- list()
ExPathwayLen <- numeric()
for (i in 1:281) {
  #t <- paste("E:/PhD/Manuscript/DATA/pathwayinformation/pathwayNo/", i, sep = '')
  t1 <- paste("E:/PhD/Manuscript/DATA/BRCA/NEW/extend/RE_N/", i, sep ='')
  t2 <- paste("E:/PhD/Manuscript/DATA/BRCA/NEW/extend/RE_T/", i, sep ='')
  text1 <- readLines(paste(t1, "/sub_node_score.noa", sep =''), encoding = "UTF-8")
  text2 <- readLines(paste(t2, "/sub_node_score.noa", sep =''), encoding = "UTF-8")
  g1<-character(0)
  g2<-character(0)
  j<-1
  for(k in 2:length(text1))
  {
    g1[j] <- strsplit(text1[k]," = ")[[1]][1]; 
    j <- j + 1;
  }
  j<-1
  for(k in 2:length(text2))
  {
    g2[j] <- strsplit(text2[k]," = ")[[1]][1]
    j <- j+1;
  }
  ExPathwayList[[i]] <- union(g1, g2)
  ExPathwayLen[i] <- length(union(g1, g2))
}

PathwayList <- list()
PathwayListLen <- numeric()
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
  PathwayListLen[i] <- length(path$V1)
}
pathwayID <- read.delim("E:/PhD/Manuscript/DATA/BRCA/Result/OldPathway/pathwayinformation/pathwayID.txt", stringsAsFactors=FALSE)
pathwaynames <- pathwayID$Symbel
Pathwaylength <- data.frame(pathwaynames, PathwayListLen, ExPathwayLen)
write.table(Pathwaylength, "Pathwaylength.txt", col.names = T,row.names = F, sep = "\t", quote = F)
