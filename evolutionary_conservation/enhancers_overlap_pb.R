#########################################################################################################################
#### Histogram with number of samples/cell types in which an interaction is observed #####
sp = "human"
setwd(paste("/home/laverre/Documents/Regulatory_Landscape/data/", sp, "/", sep=""))

simul <- read.table("Simulations/samples_simulated_all_interactions_bait_all.txt_infos", header=T)
obs <- read.table("all_interactions/all_interactions_chr_info.txt", header=T)

simul$bait_chr <- factor(simul$bait_chr, levels=levels(simul$chr))
simul$cis <- ifelse(simul$bait_chr == simul$chr, "TRUE", "FALSE")
simul <- simul[which(simul$dist < 10000000 & simul$dist > 25000),]

obs$bait_chr <- factor(obs$bait_chr, levels=levels(obs$chr))
obs$cis <- ifelse(obs$bait_chr == obs$chr, "TRUE", "FALSE")
obs <- obs[which(obs$dist < 10000000 & obs$dist > 25000),]

simul_merged <- read.table("Simulations/samples_simulated_all_interactions_bait_all_merged.txt_corrected", header=T)
obs_merged <- read.table("all_interactions/all_interactions_chr_merged.txt_cell_names_corrected2", header=T)

simul_merged$bait_chr <- factor(simul_merged$bait_chr, levels=levels(simul_merged$chr))
simul_merged$cis <- ifelse(simul_merged$bait_chr == simul_merged$chr, "TRUE", "FALSE")
simul_merged <- simul_merged[which(simul_merged$dist < 10000000 & simul_merged$dist > 25000),]

obs_merged$bait_chr <- factor(obs_merged$bait_chr, levels=levels(obs_merged$chr))
obs_merged$cis <- ifelse(obs_merged$bait_chr == obs_merged$chr, "TRUE", "FALSE")
obs_merged <- obs_merged[which(obs_merged$dist < 10000000 & obs_merged$dist > 25000),]

# Nb base pairs
obs$length <- obs$end-obs$start
obs_merged$length <- obs_merged$end - obs_merged$start
simul$length <- simul$end-simul$start
simul_merged$length <- simul_merged$end - simul_merged$start

obs$nb_type <- as.factor(obs$nb_type)
obs_merged$nb_type <- as.factor(obs_merged$nb_type)
simul$nb_type <- as.factor(simul$nb_type)
simul_merged$nb_type <- as.factor(simul_merged$nb_type)

nb_pair <- data.frame(observed = sapply(levels(obs$nb_type), function(x) 
  sum(obs[which(obs$nb_type == x & obs$cis == TRUE & obs$baited_frag == "unbaited"),]$length)
  /sum(obs[which(obs$cis == TRUE & obs$baited_frag == "unbaited"),]$length)))

nb_pair$obs_merged <- sapply(levels(obs$nb_type), function(x)
  sum(obs_merged[which(obs_merged$nb_type == x & obs_merged$cis == TRUE & obs_merged$baited_frag == "unbaited"),]$length)
  /sum(obs_merged[which(obs_merged$cis == TRUE & obs_merged$baited_frag == "unbaited"),]$length))

nb_pair$simul <- sapply(levels(obs$nb_type), function(x)
  sum(simul[which(simul$nb_type == x & simul$cis == TRUE & simul$baited_frag == "unbaited"),]$length)
  /sum(simul[which(simul$cis == TRUE & simul$baited_frag == "unbaited"),]$length))

nb_pair$simul_merged <- sapply(levels(obs$nb_type), function(x)
  sum(simul_merged[which(simul_merged$nb_type == x & simul_merged$cis == TRUE & simul_merged$baited_frag == "unbaited"),]$length)
  /sum(simul_merged[which(simul_merged$cis == TRUE & simul_merged$baited_frag == "unbaited"),]$length))

par(mfrow=c(1,2))
barplot(t(as.matrix(nb_pair[,c(1,3)])), beside=T, main="Human before merge", xlab='Cell type number', ylim=c(0,0.7),
        ylab="Base pairs proportion", col=c("forestgreen", "firebrick1"),lwd=2, cex.axis=1, cex.lab=1, cex.names=1)

legend("topright", legend=c("Observed", "Simulated"),fill=c("forestgreen", "firebrick1"), bty='n')

barplot(t(as.matrix(nb_pair[,c(2,4)])), beside=T, main="After merge", xlab='Cell number', ylim=c(0,0.7),
        ylab="Base pairs proportion", col=c("forestgreen", "firebrick1"),lwd=2, cex.axis=1, cex.lab=1, cex.names=1)

#################################################################################################################
### Overlap with enhancers ####
setwd("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/")

sp = "human"
human <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/all_interactions_chr_merged.txt_cell_names", header=T)
human <- human[,1:12]
human_simul <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/human/Simulations/simulations_human_10Mb_bin5kb_fragoverbin_chr_merged.txt", header=T, sep="\t")

human$bait_chr <- factor(human$bait_chr, levels=levels(human$chr))
human$cis <- ifelse(human$bait_chr == human$chr, "TRUE", "FALSE")
human <- human[which(human$cis == TRUE),]
human <- human[which(human$baited_frag == "unbaited"),]

human_simul$bait_chr <- factor(human_simul$bait_chr, levels=levels(human_simul$chr))
human_simul$cis <- ifelse(human_simul$bait_chr == human_simul$chr, "TRUE", "FALSE")

human$other_name <- paste(human$chr, human$start, human$end)
human_simul$other_name <- paste(human_simul$chr, human_simul$start, human_simul$end)
human$length <- human$end - human$start
human_simul$length <- human_simul$end - human_simul$start

enh_dataset <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")
for (enhancer in enh_dataset){
  overlap <- read.table(paste0("/home/laverre/Documents/Regulatory_Landscape/data/human/overlap/human_merged_overlap_", enhancer,".txt_counting_pb"), header=T)
  overlap <- overlap[which(overlap$overlap_ID != "NA"),]
  overlap$frag_name <- paste(overlap$chr, overlap$start, overlap$end)
  human[[paste0(enhancer)]] <- human$other_name %in% overlap$frag_name
  human_simul[[paste0(enhancer)]] <- human_simul$other_name %in% overlap$frag_name
  
  human[[paste0(enhancer, "_proportion")]] <- ifelse(human$other_name %in% overlap$frag_name, overlap$X.overlap, 0)
  human_simul[[paste0(enhancer, "_proportion")]] <- ifelse(human_simul$other_name %in% overlap$frag_name, overlap$X.overlap, 0)
  human[[paste0(enhancer, "_overlap")]] <- ifelse(human$other_name %in% overlap$frag_name, overlap$nb_bp_overlap, 0)
  human_simul[[paste0(enhancer, "_overlap")]] <- ifelse(human_simul$other_name %in% overlap$frag_name, overlap$nb_bp_overlap, 0)
  
}

human$dist_class <-cut(human$dist, breaks=seq(from=25000, to=5000000, by=50000), include.lowest = T)
human_simul$dist_class <- cut(human_simul$dist, breaks=seq(from=25000, to=5000000, by=50000), include.lowest = T)

human_enh_dist <- data.frame(matrix(vector(), length(levels(human$dist_class)), 1)) 
human_enh_dist_simul <- data.frame(matrix(vector(), length(levels(human$dist_class)), 1)) 

# Proportion of the sequences that is enhancer
for (enhancer in enh_dataset){
  human_enh_dist[[paste0(enhancer)]] <- sapply(levels(human$dist_class), function(x) 
    sum(human[which(human$dist_class == x),][[paste0(enhancer, "_overlap")]])
    /sum(human[which(human$dist_class == x),]$length))

  human_enh_dist_simul[[paste0(enhancer)]] <- sapply(levels(human_simul$dist_class), function(x)
    sum(human_simul[which(human_simul$dist_class == x),][[paste0(enhancer, "_overlap")]])
    /sum(human_simul[which(human_simul$dist_class == x),]$length))
}

human_enh_dist <- human_enh_dist[-1]
rownames(human_enh_dist) <- levels(human$dist_class)
human_enh_dist_simul <- human_enh_dist_simul[-1]
rownames(human_enh_dist_simul) <- levels(human_simul$dist_class)

CEX = 1.1
LWD = 1

#### Number of enhancers proportion according to genomic distance ####
plot(human_enh_dist$CAGE[0:50], type="l", col="red", lwd=LWD, ylab="Sum(overlap_enh)/Sum(length_seq)",
     xlab="Linear distance to promoters regions (pb)", xaxt = "n", cex.lab=CEX, ylim=c(0,0.2), cex.axis=CEX)

points(human_enh_dist_simul$CAGE[0:50], type="l", lty=2, col="red", lwd=0.6)
points(human_enh_dist$ENCODE[0:50], type="l", col="dodgerblue3", lwd=LWD)
points(human_enh_dist_simul$ENCODE[0:50], type="l",  lty=2, col="dodgerblue3", lwd=0.6)
points(human_enh_dist$RoadMap[0:50], type="l", col="forestgreen", lwd=LWD)
points(human_enh_dist_simul$RoadMap[0:50], type="l", lty=2,  col="forestgreen", lwd=0.6)
points(human_enh_dist$GRO_seq[0:50], type="l", col="orange", lwd=LWD)
points(human_enh_dist_simul$GRO_seq[0:50], type="l", lty=2, col="orange", lwd=0.6)

class_leg <- c("25Kb", "500Kb", "1Mb", "1.5Mb", "2Mb", "2.5Mb", "3Mb", "3.5Mb", "4Mb")
axis(1, at=seq(1,81,10), labels=F)
text(seq(1,81,10), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE, cex=CEX)
legend("topright", legend=c("RoadMap", "ENCODE", "GRO_seq","FANTOM5"), fill=c("forestgreen", "dodgerblue3", "orange","red"), bty='n', cex=CEX-0.1)


