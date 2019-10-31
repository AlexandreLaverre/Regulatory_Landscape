setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/Sequence_conservation/")

Align_score <- function(sp_origin, sp_target, data){
  dataframe <- read.table(paste("PIR_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,data,"_merged.txt", sep=""), header=T)
  
  dataframe$Allexon_ungapped <- dataframe$exclude_ungapped/dataframe$all_exclude
  dataframe[which(dataframe$Allexon_ungapped == "NaN"),]$Allexon_ungapped <- 0
  
  return(dataframe$Allexon_ungapped)
}

sp_origin = "mouse"
species <- c("opossum", "rat", "human", "rabbit", "dog", "elephant", "cow", "macaque")

####### A - Dataframe : sequence conservation target to other species ####### 
# Observed
origin_conserv_seq <- read.table(paste("PIR_cons_all_overlap_PECAN_", sp_origin,"2chicken_merged.txt", sep=""), header=T)
origin_conserv_seq$length <- origin_conserv_seq$end - origin_conserv_seq$start
origin_conserv_seq$chicken<- origin_conserv_seq$exclude_ungapped/origin_conserv_seq$all_exclude
origin_conserv_seq[which(origin_conserv_seq$chicken == "NaN"),]$chicken <- 0
origin_conserv_seq <- origin_conserv_seq[,-(10:15)]
for (sp_target in species){origin_conserv_seq[[sp_target]] <- Align_score(sp_origin, sp_target, "")}

conserv <- origin_conserv_seq[,18:26]
o <- order(apply(conserv, 2, mean), decreasing = T)
obs_conserv <- data.frame(score = apply(conserv[,o], 2, mean))
obs_conserv$int_start <- apply(conserv[,o], 2, function(x) t.test(x)[["conf.int"]][1])
obs_conserv$int_end <- apply(conserv[,o], 2, function(x) t.test(x)[["conf.int"]][2]) 

#obs_conserv <- boxplot(conserv[,o], plot=F)

# Observed with enhancers
conserv_enh <- origin_conserv_seq[which(origin_conserv_seq$CAGE_count >1),18:26]
obs_conserv_enh <- data.frame(score = apply(conserv_enh[,o], 2, mean))
obs_conserv_enh$int_start <- apply(conserv_enh[,o], 2, function(x) t.test(x)[["conf.int"]][1])
obs_conserv_enh$int_end <- apply(conserv_enh[,o], 2, function(x) t.test(x)[["conf.int"]][2]) 

# Simulated
origin_conserv_seq_simul <- read.table(paste("PIR_cons_all_overlap_PECAN_", sp_origin,"2chicken_simul_merged.txt", sep=""), header=T)
origin_conserv_seq_simul$length <- origin_conserv_seq_simul$end - origin_conserv_seq_simul$start
origin_conserv_seq_simul$chicken<- origin_conserv_seq_simul$exclude_ungapped/origin_conserv_seq_simul$all_exclude
origin_conserv_seq_simul[which(origin_conserv_seq_simul$chicken == "NaN"),]$chicken <- 0
origin_conserv_seq_simul <- origin_conserv_seq_simul[,-(10:15)]

for (sp_target in species){origin_conserv_seq_simul[[sp_target]] <- Align_score(sp_origin, sp_target, "_simul")}

conserv_simul <- origin_conserv_seq_simul[,18:26]
simul_conserv <- data.frame(score = apply(conserv_simul[,o], 2, mean))
simul_conserv$int_start <- apply(conserv_simul[,o], 2, function(x) t.test(x)[["conf.int"]][1])
simul_conserv$int_end <- apply(conserv_simul[,o], 2, function(x) t.test(x)[["conf.int"]][2]) 


####### B - Conservation ~ Distance to promoters (mouse/human) #######
sp_target = "human"
origin_conserv_seq$class <-cut(origin_conserv_seq$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
origin_conserv_seq_simul$class <-cut(origin_conserv_seq_simul$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
class_leg <- c("50Kb", "500Kb", "1Mb", "1.5Mb", "2Mb","2.5Mb")

# Observed
obs_dist <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) mean(origin_conserv_seq[which(origin_conserv_seq$class == x),][[sp_target]])))
obs_dist$int_start <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x),][[sp_target]])[["conf.int"]][1])
obs_dist$int_end <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x),][[sp_target]])[["conf.int"]][2])

# Observed with enhancers
obs_dist_enh <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) mean(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),][[sp_target]])))
obs_dist_enh$int_start <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),][[sp_target]])[["conf.int"]][1])
obs_dist_enh$int_end <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),][[sp_target]])[["conf.int"]][2])

# Simulated
simul_dist <- data.frame(inter = sapply(levels(origin_conserv_seq_simul$class), function(x) mean(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),][[sp_target]])))
simul_dist$int_start <- sapply(levels(origin_conserv_seq_simul$class), function(x) t.test(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),][[sp_target]])[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(origin_conserv_seq_simul$class), function(x) t.test(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),][[sp_target]])[["conf.int"]][2])


#######  C - Overlap with PhastCons element #######
origin_conserv_seq$phastcons_part <- origin_conserv_seq$phastcons_noexonic250 / origin_conserv_seq$length
origin_conserv_seq_simul$phastcons_part <- origin_conserv_seq_simul$phastcons_noexonic250 / origin_conserv_seq_simul$length

# Observed
obs_phastcons <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) mean(origin_conserv_seq[which(origin_conserv_seq$class == x),]$phastcons_part)))
obs_phastcons$int_start <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x),]$phastcons_part)[["conf.int"]][1])
obs_phastcons$int_end <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x),]$phastcons_part)[["conf.int"]][2])

# Observed with enhancers
obs_phastcons_enh <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) mean(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$phastcons_part)))
obs_phastcons_enh$int_start <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$phastcons_part)[["conf.int"]][1])
obs_phastcons_enh$int_end <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$phastcons_part)[["conf.int"]][2])

# Simulated
simul_phastcons <- data.frame(inter = sapply(levels(origin_conserv_seq_simul$class), function(x) mean(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$phastcons_part)))
simul_phastcons$int_start <- sapply(levels(origin_conserv_seq_simul$class), function(x) t.test(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$phastcons_part)[["conf.int"]][1])
simul_phastcons$int_end <- sapply(levels(origin_conserv_seq_simul$class), function(x) t.test(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$phastcons_part)[["conf.int"]][2])

#######  D - Overlap with PhastCons element #######
origin_conserv_seq$repeat_part <- origin_conserv_seq$repeat_pb / origin_conserv_seq$length
origin_conserv_seq_simul$repeat_part <- origin_conserv_seq_simul$repeat_pb / origin_conserv_seq_simul$length

# Observed
obs_repeat <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) mean(origin_conserv_seq[which(origin_conserv_seq$class == x),]$repeat_part)))
obs_repeat$int_start <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x),]$repeat_part)[["conf.int"]][1])
obs_repeat$int_end <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x),]$repeat_part)[["conf.int"]][2])

# Observed with enhancers
obs_repeat_enh <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) mean(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$repeat_part)))
obs_repeat_enh$int_start <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$repeat_part)[["conf.int"]][1])
obs_repeat_enh$int_end <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$repeat_part)[["conf.int"]][2])

# Simulated
simul_repeat <- data.frame(inter = sapply(levels(origin_conserv_seq_simul$class), function(x) mean(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$repeat_part)))
simul_repeat$int_start <- sapply(levels(origin_conserv_seq_simul$class), function(x) t.test(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$repeat_part)[["conf.int"]][1])
simul_repeat$int_end <- sapply(levels(origin_conserv_seq_simul$class), function(x) t.test(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$repeat_part)[["conf.int"]][2])


#######  D - Conservation between mouse and human ~ Number of cell  #######
obs_nb_cell <- boxplot(origin_conserv_seq$human~origin_conserv_seq$nb_cell, plot=F)
obs_nb_cell_enh <- boxplot(origin_conserv_seq[which(origin_conserv_seq$CAGE_count >1),]$human~origin_conserv_seq[which(origin_conserv_seq$CAGE_count >1),]$nb_cell, plot=F)

save.image("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/Fig4.Rdata")
