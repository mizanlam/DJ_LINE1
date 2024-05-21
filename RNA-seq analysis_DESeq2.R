#################################
#  RNAseq analysis with DESeq2  #
#################################

library(limma)
library(gdata)
library(gplots)
library(edgeR)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(viridis)


# Uploading raw counts table which has been uploaded to the github repository.

counts <- read_delim("DJL1_KD counts table.txt", delim = "\t",
                     col_names = c("GeneID", "CTRL_1", "CTRL_2", "CTRL_3", "CTRL_4", "DJL1_1", "DJL1_2", "DJL1_3"),
                     comment = "#") %>%
  dplyr::slice(-1) %>% mutate_at(vars(CTRL_1:DJL1_3), as.numeric) %>% filter(!grepl("-Mar", GeneID))
head(counts)

counts_df <-counts %>% select(-GeneID) %>% as.data.frame()
rownames(counts_df) <- counts$GeneID
head(counts_df)
class(counts_df)

# Normalization and heatmap using DESeq2- Read depth normalization. 

samples <- c("CTRL 1", "CTRL 2", "CTRL 3", "CTRL 4", "DJL1 2", "DJL1 3", "DJL1 4")
condition <- factor(c(1,1,1,1,2,2,2), labels = c("CTRL","DJ-L1"))
info <- as.data.frame(cbind(samples, condition))

info

matrix <- DESeqDataSetFromMatrix(countData = counts_df, colData = info, design = ~condition)
matrix <- matrix[rowSums(counts(matrix)) >= 1,]

# Normalize library
matrixR <- rlog(matrix)
exprsR <- assay(matrixR)
head(exprsR)
exprsR_tbl <- as_tibble(exprsR) %>% add_column(GeneID = rownames(exprsR), .before = 1)

pcaData <- plotPCA(matrixR, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA plot

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  scale_color_manual(values = c("grey","#838FFB"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(axis.text.x=element_text(size=12, color = "black"), axis.text.y=element_text(size=12, color = "black")) +
  theme(axis.title.x=element_text(vjust=-0.35, size=16), axis.title.y=element_text(vjust=0.35, size=16)) +
  theme(panel.background=element_rect(fill="transparent", color="black", size=1, linetype="solid"), plot.background = element_rect(fill = "transparent"), panel.grid.major=element_blank(), panel.grid.minor = element_blank()) +
  labs(x = paste0("PC1 (",percentVar[1],"% of variance)"), y = paste0("PC2 (",percentVar[2],"% of variance)"), title = "PCA") +
  ylim(c(-3,3)) + xlim(c(-15,20)) #+

