setwd("/home/laverre/Documents/Regulatory_Landscape/public_scripts/main_figures/")

human <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/all_interactions_chr_merged.txt_test_cell", header=T)
human_simul <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/human/Simulations/simulations_human_10Mb_bin5kb_fragoverbin_chr_merged.txt", header=T)
#mouse <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/mouse/all_interactions/all_interactions_chr_merged.txt_test_med", header=T)
#mouse_simul <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/mouse/Simulations/simulations_mouse_10Mb_bin5kb_fragoverbin_chr_merged.txt", header=T)

human$bait_chr <- factor(human$bait_chr, levels=levels(human$chr))
human$cis <- ifelse(human$bait_chr == human$chr, "TRUE", "FALSE")
human <- human[which(human$cis == TRUE),]
human <- human[which(human$baited_frag == "unbaited"),]

mouse$bait_chr <- factor(mouse$bait_chr, levels=levels(mouse$chr))
mouse$cis <- ifelse(mouse$bait_chr == mouse$chr, "TRUE", "FALSE")
mouse <- mouse[which(mouse$cis == TRUE),]
mouse <- mouse[which(mouse$baited_frag == "unbaited"),]

human_simul$bait_chr <- factor(human_simul$bait_chr, levels=levels(human_simul$chr))
human_simul$cis <- ifelse(human_simul$bait_chr == human_simul$chr, "TRUE", "FALSE")
mouse_simul$bait_chr <- factor(mouse_simul$bait_chr, levels=levels(mouse_simul$chr))
mouse_simul$cis <- ifelse(mouse_simul$bait_chr == mouse_simul$chr, "TRUE", "FALSE")

human$other_name <- paste(human$chr, human$start, human$end)
mouse$other_name <- paste(mouse$chr, mouse$start, mouse$end)
human_simul$other_name <- paste(human_simul$chr, human_simul$start, human_simul$end)
mouse_simul$other_name <- paste(mouse_simul$chr, mouse_simul$start, mouse_simul$end)

par(mfrow=c(2,2))
### A - % overlap with enhancers
enh_dataset <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")
proportion <- c()
conf_up <- c()
conf_low <- c()
comp_test <- c()
for (enhancer in enh_dataset){
  overlap <- read.table(paste0("/home/laverre/Documents/Regulatory_Landscape/data/human/overlap/human_merged_overlap_", enhancer,".txt"), header=T)
  overlap <- overlap[which(overlap$overlap_ID != "NA"),]
  overlap$frag_name <- paste(overlap$chr, overlap$start, overlap$end)
  human[[paste0(enhancer)]] <- human$other_name %in% overlap$frag_name
  
  x <- prop.test(x = nrow(human[which(human[[paste0(enhancer)]] == TRUE),]), n=nrow(human), p=0.5)
  proportion <- c(proportion, x$estimate)
  conf_up <- c(conf_up, x$conf.int[1])
  conf_low <- c(conf_low, x$conf.int[2])
  
  overlap <- read.table(paste0("/home/laverre/Documents/Regulatory_Landscape/data/human/overlap/human_simul_merged_overlap_", enhancer,".txt"), header=T)
  overlap <- overlap[which(overlap$overlap_ID != "NA"),]
  overlap$frag_name <- paste(overlap$chr, overlap$start, overlap$end)
  human_simul[[paste0(enhancer)]] <- human_simul$other_name %in% overlap$frag_name
  x <- prop.test(x = nrow(human[which(human_simul[[paste0(enhancer)]] == TRUE),]), n=nrow(human_simul), p=0.5)
  
  proportion <- c(proportion, x$estimate, 0)
  conf_up <- c(conf_up, x$conf.int[1], -1)
  conf_low <- c(conf_low, x$conf.int[2], -1)
  
  comp_test <- c(comp_test, prop.test(x = c(nrow(human[which(human[[paste0(enhancer)]] == TRUE),]), nrow(human[which(human_simul[[paste0(enhancer)]] == TRUE),])),
                                      n=c(nrow(human), nrow(human_simul)))$p.value)
}

# B - % overlap as a function of number of cell types
human$nb_type <- as.factor(human$nb_type)
human_enh_cell <- data.frame(matrix(vector(), length(levels(human$nb_type)), 1)) 

for (enhancer in enh_dataset){
  human_enh_cell[[paste0(enhancer)]] <- sapply(levels(human$nb_type), function(x) (nrow(human[which(human[[paste0(enhancer)]] == TRUE & human$nb_type == x),])/nrow((human[which(human$nb_type == x),]))))
  human_enh_cell[[paste0(enhancer, "_conflow")]] <- sapply(levels(human$nb_type), function(x)  (prop.test(x = nrow(human[which(human[[paste0(enhancer)]] == TRUE & human$nb_type == x),]), n=nrow(human[which(human$nb_type == x),]), p=0.5)$conf.int[1]))
  human_enh_cell[[paste0(enhancer, "_confup")]] <- sapply(levels(human$nb_type), function(x)  (prop.test(x = nrow(human[which(human[[paste0(enhancer)]] == TRUE & human$nb_type == x),]), n=nrow(human[which(human$nb_type == x),]), p=0.5)$conf.int[2]))
}

human_enh_cell <- human_enh_cell[-1]
rownames(human_enh_cell) <- levels(human$nb_type)

# C - % overlap as a function of the genomic distance
human$dist_class <-cut(human$dist, breaks=seq(from=25000, to=5000000, by=50000), include.lowest = T)
human_simul$dist_class <- cut(human_simul$dist, breaks=seq(from=25000, to=5000000, by=50000), include.lowest = T)

human_enh_dist <- data.frame(matrix(vector(), length(levels(human$dist_class)), 1)) 
human_enh_dist_simul <- data.frame(matrix(vector(), length(levels(human$dist_class)), 1)) 

for (enhancer in enh_dataset){
  human_enh_dist[[paste0(enhancer)]] <- sapply(levels(human$dist_class), function(x) (nrow(human[which(human[[paste0(enhancer)]] == TRUE & human$dist_class == x),])/nrow((human[which(human$dist_class == x),]))))
  human_enh_dist[[paste0(enhancer, "_conflow")]] <- sapply(levels(human$dist_class), function(x)  (prop.test(x = nrow(human[which(human[[paste0(enhancer)]] == TRUE & human$dist_class == x),]), n=nrow(human[which(human$dist_class == x),]), p=0.5)$conf.int[1]))
  human_enh_dist[[paste0(enhancer, "_confup")]] <- sapply(levels(human$dist_class), function(x)  (prop.test(x = nrow(human[which(human[[paste0(enhancer)]] == TRUE & human$dist_class == x),]), n=nrow(human[which(human$dist_class == x),]), p=0.5)$conf.int[2]))
  
  human_enh_dist_simul[[paste0(enhancer)]] <- sapply(levels(human_simul$dist_class), function(x) (nrow(human_simul[which(human_simul[[paste0(enhancer)]] == TRUE & human_simul$dist_class == x),])/nrow((human_simul[which(human_simul$dist_class == x),]))))
  human_enh_dist_simul[[paste0(enhancer, "_conflow")]] <- sapply(levels(human_simul$dist_class), function(x)  (prop.test(x = nrow(human_simul[which(human_simul[[paste0(enhancer)]] == TRUE & human_simul$dist_class == x),]), n=nrow(human_simul[which(human_simul$dist_class == x),]), p=0.5)$conf.int[1]))
  human_enh_dist_simul[[paste0(enhancer, "_confup")]] <- sapply(levels(human_simul$dist_class), function(x)  (prop.test(x = nrow(human_simul[which(human_simul[[paste0(enhancer)]] == TRUE & human_simul$dist_class == x),]), n=nrow(human_simul[which(human_simul$dist_class == x),]), p=0.5)$conf.int[2]))
  
}

human_enh_dist <- human_enh_dist[-1]
rownames(human_enh_dist) <- levels(human$dist_class)
human_enh_dist_simul <- human_enh_dist_simul[-1]
rownames(human_enh_dist_simul) <- levels(human_simul$dist_class)

save.image("Fig2_cell.Rdata")
