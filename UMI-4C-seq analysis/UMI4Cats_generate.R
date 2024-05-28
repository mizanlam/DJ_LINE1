##################################################
#Generating UMI4C Objects 
#This script does the following: 1) generates counts tables for two types alignment methodoliges
#2) Generates stats files and a master folder for these counts/stats files
#3) Generates the UMI4C R objects which can be loaded into "UMI4Cats_plot" to make the plots
#4) Generates bedGraph files 

#load packages
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(UMI4Cats)
library(Rbowtie2)
library(fs)
library(knitr)
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)

#define path
path <- ("C:/Users/siguy/OneDrive")

#load baits tables
setwd(file.path(path,"Documents/Ph.D/4C/Data/Samples_and_Baits/WL3420-WL3428"))
baits_genes <- read.table("baits_genes.txt", header = TRUE)
baits_13 <- read.table("baits_chr13.txt", header=TRUE)
baits_14 <- read.table("baits_chr14.txt", header=TRUE)
baits_15 <- read.table("baits_chr15.txt", header=TRUE)
baits_21 <- read.table("baits_chr21.txt", header=TRUE)
baits_22 <- read.table("baits_chr22.txt", header=TRUE)
baits_acro <- read.table("baits_acro.txt", header = TRUE)

#load samples
samples_Naive_vs_Primed <- read.table("Samples_WL3423-WL3428.txt", header = TRUE)
samples_Naive <- read.table("Samples_WL3423-WL3425.txt", header = TRUE)
samples_Primed <- read.table("Samples_WL3426-WL3428.txt", header = TRUE)

#genome indices each index built with all non acrocentric chromosomes + 1 of the acrocentrics for each object
index_path <- file.path(path, "/Documents/Ph.D/4C/Data/Reference_Genome_Index_T2T/T2T-CHM13v1.0/chm13.draft_v1.0_plusY") #used prebuilt index for v1 here
index_path_13 <- file.path(path, "/Documents/Ph.D/4C/Data/Reference_Genome_Index_T2T/T2T-CHM13v2.0_acro13/T2T-CHM13v2.0_acro13") #built my own indices with v2 for the individual ones
index_path_14 <- file.path(path, "/Documents/Ph.D/4C/Data/Reference_Genome_Index_T2T/T2T-CHM13v2.0_acro15/T2T-CHM13v2.0_acro14")
index_path_15 <- file.path(path, "/Documents/Ph.D/4C/Data/Reference_Genome_Index_T2T/T2T-CHM13v2.0_acro15/T2T-CHM13v2.0_acro15") 
index_path_21 <- file.path(path, "/Documents/Ph.D/4C/Data/Reference_Genome_Index_T2T/T2T-CHM13v2.0_acro15/T2T-CHM13v2.0_acro21")
index_path_22 <- file.path(path, "/Documents/Ph.D/4C/Data/Reference_Genome_Index_T2T/T2T-CHM13v2.0_acro15/T2T-CHM13v2.0_acro22")

#adding index paths to bait tables
baits_13$index <- index_path_13
baits_14$index <- index_path_14
baits_15$index <- index_path_15
baits_21$index <- index_path_21
baits_22$index <- index_path_22

#generating a list which contains all bait tables
acro_table <- list(baits_13, baits_14, baits_15, baits_21, baits_22)
names(acro_table) <- c("acro_13", "acro_14", "acro_15", "acro_21", "acro_22")

#digest genome
refgen <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
seqlevelsStyle(refgen) <- "UCSC"

digestGenome(
  cut_pos = 0,
  res_enz = "GATC",
  name_RE = "DpnII",
  ref_gen = refgen,
  sel_chr = NULL,
  out_path = file.path(path, "/Documents/Ph.D/4C/Data")
)

T2T_dpnii <- file.path(path, "/Documents/Ph.D/4C/Data/BSgenome.Hsapiens.NCBI.T2T-CHM13v2.0_DpnII")

#Making TxDB Object
gff <- import(file.path(path, "/Documents/Ph.D/4C/Data/T2T-CHM13v2.0_genomic_gff/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff"))
chrominfo <- getChromInfoFromNCBI("T2T-CHM13v2.0")
seqlevels(gff) <- setNames(chrominfo$SequenceName, chrominfo$RefSeqAccn) 
seqinfo(gff) <- Seqinfo(genome="T2T-CHM13v2.0")
TxDB_T2T <- makeTxDbFromGRanges(gff, taxonomyId=9606)

##################################################
#1st method to generate objects mapping to individual acrocentric chromosomes, 2nd method at the end of the script

##generating contacts files (.tsvs), loading 2*10^6 reads into memory
z=1
while (z < (length(acro_table) +1)) {
  y=1
  while (y < (nrow(samples_Naive_vs_Primed) +1)) {
    raw_dir = file.path(path, "Ph.D_repository/4C", samples_Naive_vs_Primed$Samples[y]) 
    x=1
      while (x < (length(acro_table[[z]]$Bait_ID) +1)) {
        contactsUMI4C(
        fastq_dir = raw_dir,
        wk_dir =  file.path(path,"Documents/Ph.D/4C/Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[x]), acro_table[[x]]$Bait_name[x]),
        bait_seq = acro_table[[z]]$Bait_seq[x],
        bait_pad = acro_table[[z]]$Bait_pad[x],
        res_enz = "GATC",
        cut_pos = 0,
        sel_seqname = acro_table[[z]]$Bait_chr[x],
        digested_genome = T2T_dpnii,
        bowtie_index = acro_table[[z]]$index[x],
        ref_gen = refgen,
        threads = 1,
        rm_tmp = TRUE,
        numb_reads = 2e+06,
      )
      x = x+1
    }
    y = y+1 
  }
  z = z+1
}

#2nd method aligns to the full T2T build, but allows for multiple alignments with the -k 5 flag in bowtie2

#editing alignment function, change to have "-k 5" in Rbotwtie2 to 
#add "-k 5", to line 12 after threads

trace(".singleAlignmentUMI4C", edit = TRUE, where=contactsUMI4C)

while (y < (nrow(samples_Naive_vs_Primed) +1)) {
  raw_dir = file.path(path, "Ph.D_repository/4C", samples_Naive_vs_Primed$Samples[y]) 
  x=1
  while (x < (nrow(baits_acro) +1)) {
    contactsUMI4C(
      fastq_dir = raw_dir,
      wk_dir =  file.path(path,"Documents/Ph.D/4C/Data", samples_Naive_vs_Primed$Samples[y], baits_acro$Bait_name[x]),
      bait_seq = baits_acro$Bait_seq[x],
      bait_pad = baits_acro$Bait_pad[x],
      res_enz = "GATC",
      cut_pos = 0,
      sel_seqname = baits_acro$Bait_chr[x],
      digested_genome = T2T_dpnii,
      bowtie_index = index_path,
      ref_gen = refgen, # Input bait chr to reduce running time
      threads = 1,
      rm_tmp = FALSE,
      numb_reads = 2e+06,
    )
    x = x +1
  }
  y = y +1 
}

x=1

##################################################
#Showing stat file generation, master folder, and UMI4C objects for methodology 1 

#generating stats files within the logs folder of each viewpoint
z=1
while (z < (length(acro_table)+1)) {
  y=1
  while (y < (nrow(samples_Naive_vs_Primed)+1)) {
    x = 1
    while (x < (length(acro_table[[z]]$Bait_ID) +1)) {
      setwd(file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[z]), acro_table[[z]]$Bait_name[x]))
      p = statsUMI4C(wk_dir = file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[z]), acro_table[[z]]$Bait_name[x]))
      pdf (file = paste(acro_table[[z]]$Bait_name[x], "_stats.pdf", sep = ""), width = 11, height = 8)
      print(p)
      dev.off()
      x = x+1
    }
    y = y+1
  }
  z = z+1
}

#making Master folder with renamed counts.tsv files corresponding to each viewpoint
z=1
while (z < (length(acro_table) +1)) {
  y = 1
  while (y < (nrow(samples_Naive_vs_Primed) +1)) {
    x = 1 
    dir.create(file.path(path,"Documents/Ph.D/4C/Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[z]), "Master"))
    dir.create(file.path(path,"Documents/Ph.D/4C/Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[z]), "Master", "logs"))
    dir.create(file.path(path,"Documents/Ph.D/4C/Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[z]), "Master", "count"))
    master_path <- file.path(path,"Documents/Ph.D/4C/Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[z]), "Master", "count")
    logs_path <- file.path(path,"Documents/Ph.D/4C/Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[z]), "Master", "logs")
    align <- data.frame(matrix(vector(),ncol=4))
    reads <- data.frame(matrix(vector(),ncol=4))
    stats <- data.frame(matrix(vector(), ncol=8))
    colnames(reads) <- c("sample_id", "total_reads",	"specific_reads",	"filtered_reads")
    colnames(align) <- c("sample_id",	"al_mapped",	"al_unmapped",	"al_secondary")
    colnames(stats) <- c("sample_id",	"specific_reads",	"nonspecific_reads",	"filtered_reads",	"filtout_reads",	"al_mapped",	"al_unmapped",	"umi")  
    while (x < (length(acro_table[[z]]$Bait_ID) +1)) {
      setwd(file.path(path,"Documents/Ph.D/4C/Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[z]), acro_table[[z]]$Bait_name[x], "logs"))
      align_stats <- read.delim("umi4c_alignment_stats.txt")
      read_stats <- read.delim("umi4c_stats.txt")
      umi_stats <- read.delim("stats_summary.txt")
      align[nrow(align) + 1,] <- c(paste(acro_table[[z]]$Bait_name[x], "_R1", sep = ""), align_stats$al_mapped[1], align_stats$al_unmapped[1], align_stats$al_secondary[1])
      align[nrow(align) + 1,] <- c(paste(acro_table[[z]]$Bait_name[x], "_R2", sep= ""), align_stats$al_mapped[2], align_stats$al_unmapped[2], align_stats$al_secondary[2])
      reads[nrow(reads) + 1,] <- c(acro_table[[z]]$Bait_name[x], read_stats$total_reads, read_stats$specific_reads, read_stats$filtered_reads)
      stats[nrow(stats) + 1,] <- c(acro_table[[z]]$Bait_name[x], umi_stats$specific_reads, umi_stats$nonspecific_reads, umi_stats$filtered_reads, umi_stats$filtout_reads, umi_stats$al_mapped, umi_stats$al_unmapped, umi_stats$umi) 
      source_path <- file.path(path, "Documents/Ph.D/4C/Data", samples_Naive_vs_Primed$Samples[y], names(acro_table[z]), acro_table[[z]]$Bait_name[x], "count")
      datafiles <- dir(source_path, "*.tsv.gz")
      file.copy(file.path(source_path, datafiles), master_path, overwrite = TRUE)
      file.rename(file.path(master_path, datafiles), file.path(master_path, paste(acro_table[[z]]$Bait_name[x], ".tsv.gz", sep = "")))
      x = x +1
    }
    write.table(align, file = file.path(logs_path, "umi4c_alignment_stats.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(reads, file = file.path(logs_path, "umi4c_stats.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(stats, file = file.path(logs_path, "stats_summary.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    y = y + 1
  }
  z = z +1
}

##name samples
sampleID_Naive_vs_Primed <- c("Naive_1", "Naive_2", "Naive_3", "Primed_1", "Primed_2", "Primed_3")
replicate_Naive_vs_Primed <- c(1,2,3,1,2,3)
condition_Naive_vs_Primed <-c("Naive","Naive", "Naive", "Primed", "Primed", "Primed")

##making colData files
z=1
while (z < (length(acro_table) +1)) {
  x=1
  while (x < (length(acro_table[[z]]$Bait_ID) +1)) {
    current_dir = file.path(path, "Documents", "Ph.D", "4C", "Data", "Viewpoints_WL3420-WL3428")
    viewpoint = acro_table[[z]]$Bait_name[[x]]
    file = c(file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[1], names(acro_table[z]), "Master", "count", paste(acro_table[[z]]$Bait_name[x], ".tsv.gz", sep = "")),
         file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[2], names(acro_table[z]), "Master", "count", paste(acro_table[[z]]$Bait_name[x], ".tsv.gz", sep = "")),
         file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[3], names(acro_table[z]),"Master", "count", paste(acro_table[[z]]$Bait_name[x], ".tsv.gz", sep = "")),
         file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[4], names(acro_table[z]),"Master", "count", paste(acro_table[[z]]$Bait_name[x], ".tsv.gz", sep = "")),
         file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[5], names(acro_table[z]),"Master", "count", paste(acro_table[[z]]$Bait_name[x], ".tsv.gz", sep = "")),
         file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[6], names(acro_table[z]),"Master", "count", paste(acro_table[[z]]$Bait_name[x], ".tsv.gz", sep = "")))
    colData_df = data.frame(sampleID_Naive_vs_Primed, replicate_Naive_vs_Primed, condition_Naive_vs_Primed, viewpoint, file)
    names(colData_df) = c("sampleID", "replicate", "condition", "viewpoint", "file")
    dir.create(file.path(current_dir, names(acro_table[z]),  acro_table[[z]]$Bait_name[x]))
    setwd(file.path(path,"Documents/Ph.D/4C/Data/Viewpoints_WL3420-WL3428", names(acro_table[z]), acro_table[[z]]$Bait_name[x]))
    write.table(colData_df, file = "colData_Naive_vs_primed.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    x = x +1
  }   
  z = z +1
}

#generating UMI4C objects for all acrocentric chromosomes normalized to the sample with the lowest number of UMIs
UMI_Info_chr13_1.5MB <- list()
view_dir <- file.path(path, "Documents/Ph.D/4C/Data/Viewpoints_WL3420-WL3428")
x = 1
while (x < (length(acro_table[[1]]$Bait_ID) +1)) {
  setwd(file.path(view_dir, names(acro_table[1]), acro_table[[1]]$Bait_name[x]))
  colData <- read.table("colData_Naive_vs_primed.txt", header = TRUE)
  UMI_Info_chr13_1.5MB[[x]] <- makeUMI4C(colData = colData,
                                        viewpoint_name = acro_table[[1]]$Bait_name[x],
                                        normalized = TRUE,
                                        grouping = "condition",
                                        bait_expansion = 1500000,
                                        bait_exclusion = 10000,
                                        min_win_factor = 0.02)
  x = x + 1
}

UMI_Info_chr14_1.5MB <- list()
x = 1
while (x < (length(acro_table[[2]]$Bait_ID) +1)) {
  setwd(file.path(view_dir, names(acro_table[2]), acro_table[[2]]$Bait_name[x]))
  colData <- read.table("colData_Naieve_vs_primed.txt", header = TRUE)
  UMI_Info_chr14_1.5MB[[x]] <- makeUMI4C(colData = colData,
                                        viewpoint_name = acro_table[[2]]$Bait_name[x],
                                        normalized = TRUE,
                                        grouping = "condition",
                                        bait_expansion = 1500000,
                                        bait_exclusion = 10000,
                                        min_win_factor = 0.02)
  x = x + 1
}

UMI_Info_chr15_1.5MB <- list()
x = 1
while (x < (length(acro_table[[3]]$Bait_ID) +1)) {
  setwd(file.path(view_dir, names(acro_table[3]), acro_table[[3]]$Bait_name[x]))
  colData <- read.table("colData_Naieve_vs_primed.txt", header = TRUE)
  UMI_Info_chr15_1.5MB[[x]] <- makeUMI4C(colData = colData,
                                            viewpoint_name = acro_table[[3]]$Bait_name[x],
                                            normalized = TRUE,
                                            grouping = "condition",
                                            bait_expansion = 1500000,
                                            bait_exclusion = 10000,
                                            min_win_factor = 0.02)
  x = x + 1
}

UMI_Info_chr21_1.5MB <- list()
x = 1
while (x < (length(acro_table[[4]]$Bait_ID) +1)) {
  setwd(file.path(view_dir, names(acro_table[4]), acro_table[[4]]$Bait_name[x]))
  colData <- read.table("colData_Naieve_vs_primed.txt", header = TRUE)
  UMI_Info_chr21_1.5MB[[x]] <- makeUMI4C(colData = colData,
                                        viewpoint_name = acro_table[[4]]$Bait_name[x],
                                        normalized = TRUE,
                                        grouping = "condition",
                                        bait_expansion = 1500000,
                                        bait_exclusion = 10000,
                                        min_win_factor = 0.02)
  x = x + 1
}

UMI_Info_chr22_1.5MB <- list()
x = 1
while (x < (length(acro_table[[5]]$Bait_ID) +1)) {
  setwd(file.path(view_dir, names(acro_table[5]), acro_table[[5]]$Bait_name[x]))
  colData <- read.table("colData_Naieve_vs_primed.txt", header = TRUE)
  UMI_Info_chr22_1.5MB[[x]] <- makeUMI4C(colData = colData,
                                        viewpoint_name = acro_table[[5]]$Bait_name[x],
                                        normalized = TRUE,
                                        grouping = "condition",
                                        bait_exclusion = 10000,
                                        bait_expansion = 1500000,
                                        min_win_factor = 0.02)
  x = x + 1
}

save.image("UMI4C_Objs_acro.RData") #to generate objects which can be loaded into UMI4Cats_plot.R

##################################################
#Generating bedGraph files for methodology 2 (Figure S3) 
while (x < (nrow(baits_acro)+1)) {
  y = 1
  setwd(file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[[y]], baits_acro$Bait_name[[x]], "count"))
  Naive_1 <- read.table(gzfile(list.files(pattern = "*.tsv.gz")), header = TRUE)
  y = y + 1 
  setwd(file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[[y]], baits_acro$Bait_name[[x]], "count"))
  Naive_2 <- read.table(gzfile(list.files(pattern = "*.tsv.gz")), header = TRUE)
  y = y + 1 
  setwd(file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[[y]], baits_acro$Bait_name[[x]], "count"))
  Naive_3 <- read.table(gzfile(list.files(pattern = "*.tsv.gz")), header = TRUE)
  y = y + 1 
  setwd(file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[[y]], baits_acro$Bait_name[[x]], "count"))
  Primed_1 <- read.table(gzfile(list.files(pattern = "*.tsv.gz")), header = TRUE)
  y = y + 1 
  setwd(file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[[y]], baits_acro$Bait_name[[x]], "count"))
  Primed_2 <- read.table(gzfile(list.files(pattern = "*.tsv.gz")), header = TRUE)
  y = y + 1 
  setwd(file.path(path, "Documents", "Ph.D", "4C", "Data", samples_Naive_vs_Primed$Samples[[y]], baits_acro$Bait_name[[x]], "count"))
  Primed_3 <- read.table(gzfile(list.files(pattern = "*.tsv.gz")), header = TRUE)
  Naive_avg <- Naive_1$UMIs + Naive_2$UMIs + Naive_3$UMIs
  Primed_avg <- Primed_1$UMIs + Primed_2$UMIs + Primed_3$UMIs
  Naive_bedgraph <- data.frame(Naive_1$chr_contact, Naive_1$start_contact, Naive_1$end_contact, Naive_avg)
  Primed_bedgraph <- data.frame(Primed_1$chr_contact, Primed_1$start_contact, Primed_1$end_contact, Primed_avg)
  setwd(file.path(path, "Documents", "Ph.D", "4C", "Data", "Viewpoints_WL3420-WL3428", baits_acro$Bait_name[[x]]))
  desc <- paste("track type=bedGraph name=Naive_", baits_acro$Bait_name[x], sep ="")
  writeLines(desc, paste("Naive_", baits_acro$Bait_name[[x]], ".bedGraph", sep= ""))
  write.table(Naive_bedgraph, col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t", file =paste("Naive_", baits_acro$Bait_name[[x]], ".bedGraph", sep= ""), append = T)
  desc <- paste("track type=bedGraph name=Primed_", baits_acro$Bait_name[[x]], sep ="")
  writeLines(desc, paste("Primed_", baits_acro$Bait_name[[x]], ".bedGraph", sep= ""))
  write.table(Primed_bedgraph, col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t", file =paste("Primed_", baits_acro$Bait_name[[x]], ".bedGraph", sep= ""), append = T)
  x = x +1
}
##################################################