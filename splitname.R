setwd("E:/PhD/Manuscript/DATA/BRCA/NEW/network")
library(stringr)     
first <- character()
second <- character()
for (i in 1:length(name$name)) {
  first[i] <- str_split(name$name[i],"\"")[[1]][1]
  second[i] <- str_split(name$name[i],"\"")[[1]][3]
}
newnetwork <- data.frame(first, second)
write.table(newnetwork, "newnetwork.txt", row.names = F, col.names = T, sep = "\t", quote = F )
getwd()
