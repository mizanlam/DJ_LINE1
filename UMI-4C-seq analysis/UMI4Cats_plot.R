##################################################
#Generating contact profile plots & interaction matrices requires objects made in 'UMI4Cats_generate.R'

load("UMI4C_Objs_acro.RData")

##chr15 contatcs plot primed vs naieve (Figure 2)
pdf (file = paste(acro_table[[3]]$Bait_name[1], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr15_1.5MB[[1]], dgram_plot = TRUE, xlim=c(2050000,2550000), ylim=c(0,20), rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[2], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr15_1.5MB[[2]], dgram_plot = TRUE, xlim=c(2050000,2550000), ylim=c(0,32), rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[3], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr15_1.5MB[[3]], dgram_plot = TRUE, xlim=c(2050000,2550000), ylim=c(0,10), rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[4], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr15_1.5MB[[4]], dgram_plot = TRUE, xlim=c(2050000,2550000), ylim=c(0,15), rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[5], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr15_1.5MB[[5]], dgram_plot = TRUE, xlim=c(2050000,2550000), ylim=c(0,40), rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[6], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr15_1.5MB[[6]], dgram_plot = TRUE, xlim=c(2050000,2550000), ylim=c(0,30), rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[7], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr15_1.5MB[[7]], dgram_plot = TRUE, xlim=c(2050000,2550000), ylim=c(0,42), rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

##chr13 contatcs plot primed vs naieve (Figure S2)
pdf (file = paste(acro_table[[1]]$Bait_name[1], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr13_1.5MB[[1]], dgram_plot = TRUE, xlim=c(5300000,5800000), ylim=c(0,20),  rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[2], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr13_1.5MB[[2]], dgram_plot = TRUE, xlim=c(5300000,5800000), ylim=c(0,32),  rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[3], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr13_1.5MB[[3]], dgram_plot = TRUE, xlim=c(5300000,5800000), ylim=c(0,20),  rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[4], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr13_1.5MB[[4]], dgram_plot = TRUE, xlim=c(5300000,5800000), ylim=c(0,15),  rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[5], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr13_1.5MB[[5]], dgram_plot = TRUE, xlim=c(5300000,5800000), ylim=c(0,2),  rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[6], "_naieve_vs_primed_DJ.pdf", sep = ""), width = 6, height = 6)
p <- plotUMI4C(UMI_Info_chr13_1.5MB[[6]], dgram_plot = TRUE, xlim=c(5300000,5800000), ylim=c(0,42),  rel_heights = c(0.1,0.65,0.25), colors = c("#d14b72","#909861"))
print(p)
dev.off()

##chr15 contacts plot primed vs naieve with statistical contacts for making the matrices (Figure 2 and S2)

#edit plotInteraction function, lines 14-16, colors, low = grey97, high = #92d050, na.value = grey97 
trace(plotInteractions, edit = TRUE)

pdf (file = paste(acro_table[[3]]$Bait_name[1], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 6, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr15_1.5MB[[1]], n_frags = 8)
int <- callInteractions(UMI_Info_chr15_1.5MB[[1]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr15_1.5MB[[1]], int, significant = FALSE, xlim=c(2050000,2550000), ylim=c(0,20), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[2], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags = makeWindowFragments(UMI_Info_chr15_1.5MB[[2]], n_frags = 8)
int <- callInteractions(UMI_Info_chr15_1.5MB[[2]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr15_1.5MB[[2]], int, significant = FALSE, xlim=c(2050000,2550000), ylim=c(0,32), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[3], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr15_1.5MB[[3]], n_frags = 8)
int = callInteractions(UMI_Info_chr15_1.5MB[[3]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr15_1.5MB[[3]], int, significant = FALSE, xlim=c(2050000,2550000), ylim=c(0,10), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[4], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr15_1.5MB[[4]], n_frags = 8)
int <- callInteractions(UMI_Info_chr15_1.5MB[[4]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr15_1.5MB[[4]], int, significant = FALSE, xlim=c(2050000,2550000), ylim=c(0,15), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[5], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr15_1.5MB[[5]], n_frags = 8)
int <- callInteractions(UMI_Info_chr15_1.5MB[[5]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr15_1.5MB[[5]], int, significant = FALSE, xlim=c(2050000,2550000), ylim=c(0,40), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[6], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr15_1.5MB[[6]], n_frags = 8)
int <- callInteractions(UMI_Info_chr15_1.5MB[[6]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr15_1.5MB[[6]], int, significant = FALSE, xlim=c(2050000,2550000), ylim=c(0,30), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[3]]$Bait_name[7], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr15_1.5MB[[7]], n_frags = 8)
int <- callInteractions(UMI_Info_chr15_1.5MB[[7]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr15_1.5MB[[7]], int, significant = FALSE, xlim=c(2050000,2550000), ylim=c(0,42), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()


#chr13 contacts plot primed vs naieve with statistical contacts for making the matrices
pdf (file = paste(acro_table[[1]]$Bait_name[1], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 6, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr13_1.5MB[[1]], n_frags = 8)
int <- callInteractions(UMI_Info_chr13_1.5MB[[1]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr13_1.5MB[[1]], int, significant = FALSE, xlim=c(5300000,5800000), ylim=c(0,20), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[2], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr13_1.5MB[[2]], n_frags = 8)
int <- callInteractions(UMI_Info_chr13_1.5MB[[2]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr13_1.5MB[[2]], int, significant = FALSE, xlim=c(5300000,5800000), ylim=c(0,32), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[3], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags  <-makeWindowFragments(UMI_Info_chr13_1.5MB[[3]], n_frags = 8)
int <- callInteractions(UMI_Info_chr13_1.5MB[[3]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr13_1.5MB[[3]], int, significant = FALSE, xlim=c(5300000,5800000), ylim=c(0,20), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[4], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr13_1.5MB[[4]], n_frags = 8)
int <- callInteractions(UMI_Info_chr13_1.5MB[[4]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr13_1.5MB[[4]], int, significant = FALSE, xlim=c(5300000,5800000), ylim=c(0,15), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[5], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr13_1.5MB[[5]], n_frags = 8)
int <- callInteractions(UMI_Info_chr13_1.5MB[[5]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr13_1.5MB[[5]], int, significant = FALSE, xlim=c(5300000,5800000), ylim=c(0,40), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()

pdf (file = paste(acro_table[[1]]$Bait_name[6], "_naieve_vs_primed_stat_contacts_DJ.pdf", sep = ""), width = 7, height = 6)
win_frags <- makeWindowFragments(UMI_Info_chr13_1.5MB[[6]], n_frags = 8)
int <- callInteractions(UMI_Info_chr13_1.5MB[[6]], query_regions = win_frags, padj_threshold = 0.05)
p <- plotInteractionsUMI4C(UMI_Info_chr13_1.5MB[[6]], int, significant = FALSE, xlim=c(5300000,5800000), ylim=c(0,30), rel_heights = c(0.05,0.55,0.40), colors = c("#d14b72","#909861"))
print(p)
dev.off()