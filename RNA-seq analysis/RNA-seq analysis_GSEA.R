library(fgsea)
library(tidyverse)

# Read rnk.table available on github repository 
rnk.file <- "DJL1-KD_rankedlogFC.txt"
ranks <- read.table(rnk.file,header=FALSE, colClasses = c("character", "numeric"))
colnames(ranks) <- c("ID", "logFC")
ranks <- setNames(ranks$logFC, ranks$ID)
str(ranks)

# read in gene sets of interest and save as "pathway"

List1 <- read.table("Zygote_2cell_4cell", header = F, stringsAsFactors = F)
List1 <- as.vector(List1$V1)
List2 <- read.delim("8cell_gene.gmt", header = F, stringsAsFactors = F)
List2 <- as.vector(List2$V1)
List3 <- read.table("Morulae_gene.gmt", header = F, stringsAsFactors = F)
List3 <- as.vector(List3$V1)

## continue for as many lists required (gene sets)....

pathways <- list("Zygote/2Cell/4Cell" = List1, "8Cell" = List2,"Morulae" = List3)  # Name each gene list, and list them together as pathway
str(head(pathways))

# run the fgsea function 

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=1000, nperm=1000)
head(fgseaRes)
df <- fgseaRes

# export results in table 
#write.csv(df[,-c("leadingEdge")], "embryo_gsea_results.csv")



