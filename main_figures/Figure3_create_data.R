################################################################################################################################################
library(data.table)
options(stringsAsFactors = FALSE)

source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name

ref_sp = "human" # to change human or mouse

minDistance=25e3
maxDistance=2.025e6

enhancers = c("FANTOM5", "ENCODE")
target_sp = "human"
closest_sp="rat"

if(ref_sp == "human"){
  enhancers <- c(enhancers, "RoadmapEpigenomics", "FOCS_GRO_seq")
  target_sp = "mouse"
  closest_sp="macaque"
  }

path_evol <- paste(pathFinalData, "SupplementaryDataset7/", ref_sp, "/", sep="")
path_annot <- paste(pathFinalData, "SupplementaryDataset4/", ref_sp, "/", sep="")


################################################################################################################################################
################################## Fig3-B Alignment score of contacted sequences between all species ###########################################

#### Contacted restriction fragments
obs <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp, "/statistics_contacted_sequence_original.txt", sep=""), header=T)
simul <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp,"/statistics_contacted_sequence_simulated.txt", sep=""), header=T)

obs$ID <-  do.call(paste,c(obs[c("chr","start","end")],sep=":"))
simul$ID <-  do.call(paste,c(simul[c("chr","start","end")],sep=":"))

# Filters
obs <- obs[which(obs$baited == "unbaited"),]
simul <- simul[which(simul$baited == "unbaited"),]
obs <- obs[which(obs$BLAT_match < 2),]
simul <- simul[which(simul$BLAT_match < 2),]

#### Alignment score vs all species 
align <- read.table(paste(path_evol,"sequence_conservation/restriction_fragments/Alignments_stats_all_species_nonexonic_ungapped.txt", sep="/"), header=T)
species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")

if (ref_sp=="human"){ID="ID.human"}else{ID="ID.mouse"}
align_obs <- align[which(align[,1] %in% obs$ID), c(ID, species)]
align_simul <- align[which(align[,1] %in% simul$ID), c(ID, species) ]
align_obs_enh <- align[which(align[,1] %in% obs[which(obs$ENCODE_bp > 0),]$ID), c(ID, species)] # fragment that contain at least 1 ENCODE enhancer


################################################################################################################################################
################################################# Fig3-C & G - Enhancers Alignment  ############################################################ 
list_conserv_enh <- list()
for (enh in enhancers){
  align <- read.table(paste(path_evol,"sequence_conservation/enhancers/", enh, "Alignments_stats_all_species_nonexonic_ungapped.txt", sep="/"), header=T)
  obs_stats <- read.table(paste(path_annot, enh, "/statistics_contacted_enhancers_original.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_annot, enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), header=T)
  
  obs_stats$enh <-  do.call(paste,c(obs_stats[c("chr","start","end")],sep=":"))
  simul_stats$enh <-  do.call(paste,c(simul_stats[c("chr","start","end")],sep=":"))
  
  ##### Filters
  species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
  align <- align[, c("enh", species)]
  
  obs_stats$dist_class <- cut(obs_stats$med_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  simul_stats$dist_class <-cut(simul_stats$med_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  
  # Select unduplicated and with repeat_part < 20%
  obs_stats <- obs_stats[which(obs_stats$BLAT_match < 2),] # & (obs_stats$repeat_pb/obs_stats$length) < 0.2
  simul_stats <- simul_stats[which(simul_stats$BLAT_match < 2),] # & (simul_stats$repeat_pb/simul_stats$length) < 0.2
  
  #### Fig 3-C - Distribution of alignment score of ENCODE enh according to all species
  if (enh == "ENCODE"){
    align_enhancers_obs <- align[which(align$enh %in% obs_stats$enh),]
    align_enhancers_simul <- align[which(align$enh %in% simul_stats$enh),]
  }
  
  #### Fig 3-G - Alignment score according to distance from promoters
  align_enhancers_dist <- data.frame(matrix(vector(), length(levels(obs_stats$dist_class))))
  
  align_enhancers_dist[[paste0(enh)]] <- sapply(levels(obs_stats$dist_class), function(x) 
    mean(align[which(align$enh %in% obs_stats[which(obs_stats$dist_class == x),]$enh),][[target_sp]]))
  
  align_enhancers_dist[[paste0(enh, "_start")]]  <- sapply(levels(obs_stats$dist_class), function(x) 
    t.test(align[which(align$enh %in% obs_stats[which(obs_stats$dist_class == x),]$enh),][[target_sp]])[["conf.int"]][1])
  
  align_enhancers_dist[[paste0(enh, "_end")]]  <- sapply(levels(obs_stats$dist_class), function(x) 
    t.test(align[which(align$enh %in% obs_stats[which(obs_stats$dist_class == x),]$enh),][[target_sp]])[["conf.int"]][2])
  
  list_conserv_enh <- c(list_conserv_enh, align_enhancers_dist)
}

################################################################################################################################################
################################# Fig3-D - Conservation of contacted sequences between human and mouse #########################################

obs <- obs[order(obs$ID), ]
simul <- simul[order(simul$ID), ]

obs$align_target <- align_obs[order(align_obs[,1]), target_sp]
simul$align_target <- align_simul[order(align_simul[,1]), target_sp]

obs_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$align_target)))
obs_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][1])
obs_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][2])

obs_enh_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)))
obs_enh_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][1])
obs_enh_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$align_target)))
simul_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][2])

################################################################################################################################################
################################## Fig3-E - Conservation of contacted sequences between closest species ########################################

obs$align_target <- align_obs[order(align_obs[,1]), closest_sp]
simul$align_target <- align_simul[order(align_simul[,1]), closest_sp]

obs_dist_mac <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$align_target)))
obs_dist_mac$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][1])
obs_dist_mac$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][2])

obs_enh_dist_mac <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)))
obs_enh_dist_mac$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][1])
obs_enh_dist_mac$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][2])

simul_dist_mac <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$align_target)))
simul_dist_mac$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][1])
simul_dist_mac$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][2])

################################################################################################################################################
################################## Fig3-F Proportion of repeated nucleotide according to distance ###########################################

# Proportion of non exonic repated nucleotide
obs$repeat_part <- obs$repet_noexon_bp/obs$length
simul$repeat_part <- simul$repet_noexon_bp/simul$length

# Repeat prop according to distance
obs$dist_class <- cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
simul$dist_class <-cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)

obs_repet_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$repeat_part)*100))
obs_repet_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
obs_repet_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)

obs_enh_repet_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)*100))
obs_enh_repet_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)[["conf.int"]][1]*100)
obs_enh_repet_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)[["conf.int"]][2]*100)

simul_repet_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$repeat_part)*100))
simul_repet_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
simul_repet_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)


################################################################################################################################################
########################################### Supp.Figure X -  Exonic proportion according to distance ##########################################

obs$exonic_part <- obs$all_exon_bp/obs$length
simul$exonic_part <- simul$all_exon_bp/simul$length

obs_exon_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$exonic_part)*100))
obs_exon_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$exonic_part)[["conf.int"]][1]*100)
obs_exon_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$exonic_part)[["conf.int"]][2]*100)

obs_enh_exon_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$exonic_part)*100))
obs_enh_exon_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$exonic_part)[["conf.int"]][1]*100)
obs_enh_exon_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$exonic_part)[["conf.int"]][2]*100)

simul_exon_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$exonic_part)*100))
simul_exon_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$exonic_part)[["conf.int"]][1]*100)
simul_exon_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$exonic_part)[["conf.int"]][2]*100)


################################################################################################################################################
# Output
save(align_obs, align_obs_enh, align_simul,
     obs_repet_dist, obs_enh_repet_dist, simul_repet_dist,
     obs_exon_dist, obs_enh_exon_dist, simul_exon_dist,
     obs_dist, obs_enh_dist, simul_dist,
     obs_dist_mac,obs_enh_dist_mac,simul_dist_mac,
     align_enhancers_obs, align_enhancers_simul, list_conserv_enh, 
     file = paste(pathFigures, "/Fig3_", ref_sp, ".Rdata", sep=""))


