################################
#  RNAseq analysis with EdgeR  #
################################

library(tidyverse)
library(edgeR)
library(ggplot2)
library(ggrepel)

# 1. Formatting counts tables as needed for downstream analysis

# Uploading raw counts table which has been uploaded to the github repository.

Exp <- read.table("DJL1_KD counts table.txt", row.names = 1)

d0 <- DGEList(Exp)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
d <- calcNormFactors(d)
N <- colSums(Exp)
dim(d) #number of genes left

# Adding information about samples for differential analysis

ttmt <- c(rep("CTRL",4), rep("DJ_L1",3 ))
group <- interaction(ttmt) 

# Building the model for limma-voom Read depth normaliztion  

mm <- model.matrix(~0 + group)
y <- voom(d, mm, lib.size = N, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))


# 2. Differential gene analysis with edgeR

# Differential analysis for DJ-L1 KD cells

contr <- makeContrasts(groupDJ_L1 - groupCTRL, levels = colnames(coef(fit))) 
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.01 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
DJ_L1top.table <- top.table
DJ_L1top.table$Gene <- row.names(DJ_L1top.table)


# Exporting the results
#write_csv(DJ_L1top.table, file = "export_folder/DJ_L1top.table.csv")

# 3. Volcano plot of differential gene expression

# set Categorical Variables 

DJ_L1top.table$Cat <- 0
DJ_L1top.table[DJ_L1top.table$logFC > 0.6 & DJ_L1top.table$adj.P.Val < 0.01, "Cat" ] <- 1
DJ_L1top.table[DJ_L1top.table$logFC < -0.6 & DJ_L1top.table$adj.P.Val < 0.01, "Cat" ] <- 2

# Volcano plot 

ggplot(DJ_L1top.table, aes(x=logFC, y=-log10(P.Value), col=as.factor(Cat))) +
  geom_point() +
  scale_color_manual(values = c("grey80", "#838FFB","navy")) +
  geom_hline(yintercept = 2.33,
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(1.5), -log2(1.5)),
             linetype = "dashed") +
  geom_point(data = subset(DJ_L1top.table, DJ_L1top.table$Gene %in% c("DDIT3", "CCNA1", "TRIM43", "SUSD2","H3Y1","ODC1","KDM4E","ZSCAN5B", "TPRX1","MBD3L2", "TRIM43B","RFPL2","ZSCAN4","SLC34A2","KLF17")), color = "#702dc7") + #highlightselected 8C genes 
  geom_text_repel (data = subset(DJ_L1top.table, DJ_L1top.table$Gene %in% c("DDIT3", "CCNA1","TRIM43","SUSD2","H3Y1","ODC1","KDM4E","ZSCAN5B", "TPRX1","MBD3L2", "TRIM43B","RFPL2","ZSCAN4","SLC34A2","KLF17")), 
                   color = "black",
                   segment.color = "grey50",
                   segment.linetype = 3,
                   nudge_x = 1,
                   nudge_y = 2, aes(label=Gene)) +  #to label selected 8C genes on plot
  labs(x = "Log2FC", y = "-log10(P.Value)") +
  theme(axis.text = element_text(size = 20)) +
  theme_classic()



