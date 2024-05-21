# Intro 

## ATAC-seq in Control vs DJ-L1 KD with hg38 alignment and analysis of open chromatin peaks across genome and genes
  # hg38 indexes are retrieved from UCSC   
  # Deduplicate and Index with SAMtools 
  # Call peaks with MACS3
  # Summarize reports with MultiQC. 

###### alignment_hg38.sh ######



#!/bin/bash
#SBATCH --job-name=alignment_hg38
#SBATCH -N 1
#SBATCH -n 6


# Define path to genome indices
genome_index=/home/ltri/santoslab/share/genome_files/bowtie2_indexes/hg38/hg38


for file in *_L002_R1_001.fastq.gz; do
sample=${file%_L002_R1_001.fastq.gz}


echo $sample

bowtie2 -p 6 --local -X 1000 -x $genome_index -1 $file -2 ${sample}_L002_R2_001.fastq.gz -S ${sample}.sam 2>${sample}.align.stats

echo "wc -l aligned"
wc -l ${sample}.sam

# Count reads assigned to chrM, Un, random, then remove
echo "wc -l chrM"
grep chrM ${sample}.sam | wc -l
echo "wc -l chrUn"
grep chrUn ${sample}.sam | wc -l
echo "wc -l random"
grep random ${sample}.sam | wc -l

sed '/chrM/d;/random/d;/chrUn/d' < ${sample}.sam > ${sample}removedchrs.sam

# Convert to BAM, sort
samtools view -bu ${sample}removedchrs.sam | samtools sort -o sorted_${sample}.bam

# Remove duplicates - don't need any of these files later, so name intermediate files redundantly
samtools collate -o namecollate.bam sorted_${sample}.bam
samtools fixmate -m namecollate.bam fixmate.bam
samtools sort -o positionsort.bam fixmate.bam
samtools markdup -rs positionsort.bam filtered_${sample}.bam
#markdup flag the duplicates and -r argument remove them

# Index bam
samtools index filtered_${sample}.bam

# For completeness sake, get a line count of the number of reads passing all filters, then remove the file
samtools view filtered_${sample}.bam -o filtered_${sample}.sam
echo "wc -l filtered"
wc -l filtered_${sample}.sam
rm filtered_${sample}.sam

# Remove intermediate sam files	
rm ${sample}removedchrs.sam 	
rm namecollate.bam
rm fixmate.bam
rm positionsort.bam

done

# Compile alignment statistics with multiqc
multiqc . --ignore "FastQC/"
exit 0



##### peakcalling.sh #####


#!/bin/bash
#SBATCH --job-name=peakcalling

for file in filtered_*.bam; do

# Make life easier by making the file names simpler and into an extension-less variable
sample=${file:9}
sample2=${sample%.bam}

macs3 callpeak --nomodel -B -f BAMPE -t $file -n $file --nolambda --keep-dup all --gsize hs
done

exit 0


# blacklisting peaks

#retrieved blacklist bed file from https://dozmorovlab.github.io/excluderanges/.

#!/bin/bash
#SBATCH --job-name=blacklisted

# Define path to blacklist
blacklist=/home/ltri/santoslab/share/genome_files/hg38.Kundaje.GRCh38_unified_blacklist.bed


# Remove peak file header
# Rearrange columns to bed format (chr start end name score strand), spoof strand column

# ! Check number of lines to remove from header and adjust code accordingly !
for file in *_peaks.xls;
do
awk 'NR > 26 { print }' < $file > temp
awk -v OFS='\t' {'print $1,$2,$3,$10,$6,"+"'} temp > ${file%.xls}.bed
done

rm temp

# Count number of peaks called in each sample
wc -l *peaks.bed > peak_counts

# Remove the blacklisted regions from peaks and quantify number of peaks lost
# Use bedtools to intersect peaks with blacklisted regions and exclude

for file in *peaks.bed; do
bedtools intersect -v -wa -a $file -b $blacklist > ${file%.bed}_blacklisted.bed
done

wc -l *blacklisted* > blacklisted_counts

exit 0


##### Differential peak analysis  #####

library(DiffBind) #v13.12.0

SampleID <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4", "DJL1_1", "DJL1_2", "DJL1_3","DJL1_4")
Condition <- c(rep("Control", 4), rep("DJL1", 4))
Replicate <- c(rep(seq(1:4), 1), rep(seq(1:4), 1))
bamReads <- c("filtered_Control_1_S5.bam", "filtered_Control_2_S6.bam", "filtered_Control_3_S7.bam", "filtered_Control_4_S8.bam", "filtered_DJL1_1_S1.bam", "filtered_DJL1_2_S2.bam", "filtered_DJL1_3_S3.bam", "filtered_DJL1_4_S4.bam")
Peaks <- c("Control_1_S5_peaks_blacklisted.bed", "Control_2_S6_peaks_blacklisted.bed", "Control_3_S7_peaks_blacklisted.bed", "Control_4_S8_peaks_blacklisted.bed", "DJL1_1_S1_peaks_blacklisted.bed", "DJL1_2_S2_peaks_blacklisted.bed", "DJL1_3_S3_peaks_blacklisted.bed", "DJL1_4_S4_peaks_blacklisted.bed")
PeakCaller <- c(rep("bed", 8))
sampleData_comparison <- data.frame(cbind(SampleID, Condition, Replicate, bamReads, Peaks, PeakCaller))

pdf("comparison.pdf")


comparison <- dba(sampleSheet=sampleData_comparison, minOverlap = 2)

dba.plotHeatmap(comparison, main = "No read counts - just peak presence/absence")
dba.plotPCA(comparison, DBA_CONDITION, label = DBA_ID)


comparison_counts <- dba.count(comparison, minOverlap = 2, summits = 100, bRemoveDuplicates = TRUE, score=DBA_SCORE_READS)
save.image("Export.RData")

dba.plotHeatmap(comparison_counts, score = DBA_SCORE_READS)
dba.plotPCA(comparison_counts, score = DBA_SCORE_READS, DBA_CONDITION, label = DBA_ID)


comparison_normalized <- dba.normalize(comparison_counts, library = "full", normalize = "lib",  method=DBA_DESEQ2)
dba.plotHeatmap(comparison_normalized, score = DBA_SCORE_NORMALIZED)
dba.plotPCA(comparison_normalized, score = DBA_SCORE_NORMALIZED, DBA_CONDITION, label = DBA_ID)
```

# Plot Venn diagram for Control and DJ-L1 KD separately
dba.plotVenn(comparison, comparison$masks$Control)
dba.plotVenn(comparison, comparison$masks$DJL1)

venn <- dba.peakset(comparison, consensus=DBA_CONDITION, minOverlap=2)
dba.plotVenn(venn, venn$masks$Consensus)


comparison_contrasted <- dba.contrast(comparison_normalized, contrast=c("Condition","DJL1","Control"))

comparison_analyzed <- dba.analyze(comparison_contrasted, method = DBA_DESEQ2)
dba.analyze(comparison_analyzed, bRetrieveAnalysis=TRUE)
save.image("Analysis.RData")


dba.plotMA(comparison_counts, contrast=list(DJL1=comparison_counts$masks$DJL1), bNormalized=FALSE,
           sub="Non-Normalized")

dba.plotMA(comparison_analyzed,method=DBA_DESEQ2)
dba.plotMA(comparison_analyzed,method=DBA_DESEQ2, bSignificant = FALSE)
dba.plotVolcano(comparison_analyzed, method=DBA_DESEQ2)

comparison.report <- dba.report(comparison_analyzed, th = 1, fold = 0, method = DBA_DESEQ2, bCounts = T)
comparison.df <- as.data.frame(comparison.report)

write.table(comparison.df, "comparison_DAR_table.txt", quote = F, row.names = F)


with(comparison.df, plot(Fold,-log10(FDR),
                         pch=20, main="Differential peak accessibility between Control and DJL1",
                         col = ifelse(FDR<0.05 & abs(Fold)>1, "996be3", "#C1C1C1"),
                         cex = 0.75))
abline(v=c(-1,1), lty = "dotted")
abline(h=-log10(0.05), lty = "dotted")

dev.off()


comparison.df$name <- paste0(comparison.df$seqnames, ":", comparison.df$start, "-", comparison.df$end)
# Find peaks sig up and down in DJ-L1 KD
upDJ <- subset(comparison.df, Fold > 0.75 & FDR <  0.05)
downDJ <- subset(comparison.df, Fold < -0.75 & FDR <  0.05)
dim(upDJ)

dim(downDJ)


# For viewing in IGV, sort by FC within the sig peaks to find the highest impacted regions
sorted_upDJ <- upDJ[order(-upDJ$Fold),]
sorted_downDJ <- downDJ[order(-downDJ$Fold),]

# Write out and submit to GREAT to get GO terms of nearby genes
# write.table(sorted_upDJ, "DAR_upDJ.txt", quote = F, row.names = F)
# write.table(sorted_downDJ, "DAR_downDJ.txt", quote = F, row.names = F)


######### Figures ##########

## setting categories ##

comparison.df2$diffexpressed <- "NO" #labelling all genes as NO
comparison.df2$diffexpressed[comparison.df2$FDR<0.05 & comparison.df2$Fold > 0.58 ] <- "UP"
comparison.df2$diffexpressed[comparison.df2$FDR<0.05 & comparison.df2$Fold < -0.58 ] <- "DOWN"

## Volcano Plot ##

library(ggplot2)
ggplot(data = comparison.df2 , aes(x= Fold , y= -log10(FDR), col = diffexpressed))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = c("purple", "#C1C1C1","#ff3399"))+
  geom_vline(xintercept=c(-0.58,0.58), lty="dotted") +
  geom_hline(yintercept=-log10(0.05), lty="dotted")


with(comparison.df2, plot(Fold,-log10(FDR),
                          pch=20, main="Differential peak accessibility between Control and DJL1",
                          col = ifelse(FDR<0.05 & Fold>0.58, "magenta", "#C1C1C1"),
                          cex = 0.75))
abline(v=c(-0.58,0.58), lty = "dotted")
abline(h=-log10(0.05), lty = "dotted")


## PCA Plot ##

dba.plotPCA(comparison_normalized, score = DBA_SCORE_NORMALIZED, DBA_CONDITION)


dev.off()




# Motif analysis 
  # to find transcription factor signatures associated with genomic regions preferentially accessible or inaccessible in DJ-L1 KD

## Running gkmSVM

Step 1: Identify peaks in the positive and negative training sets
Positive: Differentially OPEN regions (open in DJ-L1 KD) "BED_Up.txt"
Negative: Differentially CLOSED regions (open in Control/closed in DJ-L1 KD) "BED_Down.txt"
Input files:
  1) Positive peak.bed
2) Negative peak.bed
3) nr10mers.fa (Table with all the possible 10mer)
4) fasta files of both input peak files


## Make fasta files
```{r, eval = F}
#Make BED files starting from DAR files Up and down in DJ 

BED_Up <- sorted_upDJ_log15[, c(1,2,3)]
write.table(BED_Up, "BED_UP.bed", col.names  = FALSE, quote = FALSE, sep = "\t", row.names= FALSE)
BED_Down <- sorted_downDJ_log15[, c(1,2,3)]
write.table(BED_Down, "BED_down.bed", col.names  = FALSE, quote = FALSE, sep = "\t", row.names= FALSE)


#Open the cluster, create a folder called GKMSVM and move here the 2 new files ("BED_UP.bed" and "BED_down.bed") + the file nr10mers (given bt Sarah, it is a list of all possible combinations of 10mer)
mkdir GKMSVM
# Make fasta files

bedtools getfasta -fo Up.fa -fi /home/ltri/santoslab/share/genome_files/genome_fastas/hg38.fa -bed BED_UP.bed
bedtools getfasta -fo Down.fa -fi /home/ltri/santoslab/share/genome_files/genome_fastas/hg38.fa -bed BED_down.bed
#Check it worked (it should have the double of the number of lines of the BED file  and start with the fasta header)
wc -l Up.fa
wc -l Down.fa
less Up.fa
less Down.fa
```
## Train gkmSVM (support vector machine) in assigning enrichment scores to all the possible 10mers
```{r, eval = F}
nano gkmSVM_enriched.R
#write inside the nano
library(gkmSVM)
library(BSgenome.Hsapiens.UCSC.hg38.masked)

#take every possible 10mers and see how many times it is present in my 2 files and print it in a .out file

gkmsvm_kernel("Up.fa", "Down.fa", "Enriched_kernel.out")
#training the machine learning algorithm: take out 10% data and use them to determine how every 10mere is predictive in assigning a peak to control or ink or both conditions
#svmfnxprfx = prefix, all output file will start with this word

gkmsvm_trainCV('Enriched_kernel.out','Up.fa','Down.fa',svmfnprfx='Enriched_', outputCVpredfn='Enriched_cvpred_.out', outputROCfn='Enriched_roc.out', outputPDFfn = "Enriched_ROC.pdf")
#Assing final score of each 10mers as a result of the training (1 score for each 10mer, called weight of a 10mer)

gkmsvm_classify('nr10mers.fa',svmfnprfx='Enriched_', 'Up_vs_Down_weights.out')

nano gkmSVM_launcher.sh


#!/bin/bash
#SBATCH --job-name=gkmSVM_launcher
Rscript --save gkmSVM_enriched.R

exit 0

# Ctrl + O, save the gkmSVM_launcher.sh file
sbatch gkmSVM_launcher.sh

sort -grk2,2 Up_vs_Down_weights.out > sorted_Up_vs_Down_weights.out

R

weights <- read.table("sorted_Up_vs_Down_weights.out")
temp_wg <- as.data.frame(weights$V2)

#temporary table just for graph purpose

colnames(temp_wg) <- "scores"
temp_wg$vocabulary <- "vocabulary"

library(ggplot2)
pdf("distribution_of_scores.pdf")
ggplot(temp_wg, aes(x=vocabulary, y=scores)) +
  geom_violin() +
  geom_hline(yintercept=0.5, color='yellow') +
  geom_hline(yintercept=-0.5, color='yellow')+
  geom_hline(yintercept=1, color='orange') +
  geom_hline(yintercept=-1, color='orange')+
  geom_hline(yintercept=1.5, color='red') +
  geom_hline(yintercept=-1.5, color='red')
dev.off()

# Let's start at 0.3 as the cut-off
# Extract just those kmers
# Positive values will be enriched in DJ-L1 KD, negative values enriched in Control

weights_gainedDJ <- weights[weights$V2 > 0.3, ]
weights_lostDJ <- weights[weights$V2 < -0.3, ]

# Write out the kmers
write.table(weights_gainedDJ$V1, "kmers_gainedDJ", col.names = F, row.names = F, quote = F)
write.table(weights_lostDJ$V1, "kmers_lostDJ", col.names = F, row.names = F, quote = F)
q()

#Exit R
# Starcode to cluster kmers
#Play around with -d argument to determine the similarity inside the cluster. A lower argument value will correspond to more stringency to accept kmer in the cluster, reducing the number of kmer per cluster. We would expect a value lower than 1000 (not to high numbers). Try -d 4 (A little less restrictive than usual, but I have very few kmers for the kmers_Ink table)
#Have a look at my initial number of kmers
wc -l kmers_gainedDJ

wc -l kmers_lostDJ

starcode -d 4 -s -i kmers_gainedDJ --seq-id -o kmers_gainedDJ.clusters
starcode -d 4 -s -i kmers_lostDJ --seq-id -o kmers_lostDJ.clusters
#Have a look at the number of kmers clusters
wc -l kmers_gainedDJ.clusters

wc -l kmers_lostDJ.clusters

# From clusters to fastas for clustalomega
R
library(reshape2)
# Preparing all kmers, not clustered
kmers_gainedDJ <- read.table("kmers_gainedDJ")
kmers_lostDJ <- read.table("kmers_lostDJ")

#Attribute a number for each kmers (from 1 to N Kmers)
kmers_gainedDJ$value <- 1:nrow(kmers_gainedDJ)
kmers_lostDJ$value <- 1:nrow(kmers_lostDJ)

# Preparing the clustered kmers
## Read in
clustered_gainedDJ <- read.table(text = gsub(",", "\t", readLines("kmers_gainedDJ.clusters")), fill=TRUE)
clustered_lostDJ <- read.table(text = gsub(",", "\t", readLines("kmers_lostDJ.clusters")), fill=TRUE)

## Remove the column that just has the count of kmers in the cluster
clustered_gainedDJ <- clustered_gainedDJ[, -2]
clustered_lostDJ <- clustered_lostDJ[, -2]

## Preparing for merge
melt_clustered_gainedDJ <- melt(clustered_gainedDJ)
melt_clustered_lostDJ <- melt(clustered_lostDJ)

# Remove rows with NAs (that we inserted when making the cluster table with the argument fill=TRUE)
melt_clustered_gainedDJ <- melt_clustered_gainedDJ[complete.cases(melt_clustered_gainedDJ), ]
melt_clustered_lostDJ <- melt_clustered_lostDJ[complete.cases(melt_clustered_lostDJ), ]

# Merge the two together
merged_gainedDJ <- merge(melt_clustered_gainedDJ, kmers_gainedDJ, by="value")
merged_lostDJ <- merge(melt_clustered_lostDJ, kmers_lostDJ, by="value")
merged_gainedDJ <- merged_gainedDJ[, c(2,4)]
merged_lostDJ <- merged_lostDJ[, c(2,4)]


# Divide up into one file per centroid
## Add column names to help with split
colnames(merged_gainedDJ) <- c("centroid", "kmers")
colnames(merged_lostDJ) <- c("centroid", "kmers")

# Split by centroid
split_merged_gainedDJ <- split(merged_gainedDJ, merged_gainedDJ$centroid)
split_merged_lostDJ <- split(merged_lostDJ, merged_lostDJ$centroid)

for(centroid in names(split_merged_gainedDJ))
  write.table(split_merged_gainedDJ[[centroid]], file = paste0(centroid, "_gainedDJ.centroid"), sep="\t", row.names = F, col.names = F, quote = F)
for(centroid in names(split_merged_lostDJ))
  write.table(split_merged_lostDJ[[centroid]], file = paste0(centroid, "_lostDJ.centroid"), sep="\t", row.names = F, col.names = F, quote = F)
q()
# Turn into fastas (remove the first column, append each kmer with a ">")
for file in *.centroid; do
awk '{print $2}' $file > temp
awk '1;!(NR%1){print ">"NR;}' temp > temp1
sed -i '1i>0' temp1
sed -i '$ d' temp1
mv temp1 $file.fa
done
rm temp


#Clustalo perform local alignment on the kmers composing the cluster
mkdir clustalomega
for file in *.centroid.fa; do
clustalo -i $file -o ${file%.centroid.fa}.clustal.fa
mv ${file%.centroid.fa}.clustal.fa clustalomega
done
#Meme identify the consenus sequence for each cluster
# Loop through files to change dashes to X in prep for meme (use to make PWM logos)
mkdir meme_logos
for file in clustalomega/*.clustal.fa; do
tr - X < $file > ${file%.clustal.fa}.clustal.X
mv ${file%.clustal.fa}.clustal.X meme_logos
done

# Perform meme on all files!
for file in meme_logos/*.clustal.X; do
meme $file -maxsize 1000000 -o ${file%.clustal.X} -dna -mod zoops -nmotifs 1 -evt inf -wnsites 0.8 -minw 8 -maxw 12 -wg 11 -ws 1 -maxiter 50 -distance 0.001 -prior dirichlet -b 0.01 -spmap uni -spfuzz 0.5

done
