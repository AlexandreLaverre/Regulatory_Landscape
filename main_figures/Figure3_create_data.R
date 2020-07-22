library(data.table)
options(stringsAsFactors = FALSE)

ref_sp = "human"
target_sp = "mouse"

path <- "/home/laverre/Data/Regulatory_landscape/result/"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, "restriction_fragments", sep="/")

####################################### Restriction fragments datas  ####################################### 
path_data <- "/home/laverre/Manuscript/"
obs <- read.table(paste(path_data, "SupplementaryDataset5/", ref_sp, "/statistics_contacted_sequence_original.txt", sep=""), header=T)
simul <- read.table(paste(path_data, "SupplementaryDataset5/", ref_sp,"/statistics_contacted_sequence_simulated.txt", sep=""), header=T)

obs <- obs[which(obs$baited == "unbaited"),]
simul <- simul[which(simul$baited == "unbaited"),]
obs <- obs[which(obs$BLAT_match < 2),]
simul <- simul[which(simul$BLAT_match < 2),]

obs$ID <-  do.call(paste,c(obs[c("chr","start","end")],sep=":"))
simul$ID <-  do.call(paste,c(simul[c("chr","start","end")],sep=":"))


################# Alignment score vs all species #################
align <- read.table(paste(path_evol,"sequence_conservation/restriction_fragments/Alignments_stats_all_species_total_ungapped.txt", sep="/"), header=T)

species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")

align_obs <- align[which(align$enh %in% obs$ID), c("enh", species)]
align_obs_enh <- align[which(align$enh %in% obs[which(obs$ENCODE_bp > 0),]$ID), c("enh", species)]
align_simul <- align[which(align$enh %in% simul$ID), c("enh", species) ]


################# Repeat proportion  #################
obs$repeat_part <- obs$repet_noexon_bp/obs$length
simul$repeat_part <- simul$repet_noexon_bp/simul$length

# Repeat prop according to distance
obs$dist_class <- cut(obs$median_dist, breaks=seq(from=25000, to=2025000, by=50000), include.lowest = T)
simul$dist_class <-cut(simul$median_dist, breaks=seq(from=25000, to=2025000, by=50000), include.lowest = T)

obs_repet_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$repeat_part)*100))
obs_repet_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
obs_repet_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)

obs_enh_repet_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)*100))
obs_enh_repet_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)[["conf.int"]][1]*100)
obs_enh_repet_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)[["conf.int"]][2]*100)

simul_repet_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$repeat_part)*100))
simul_repet_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
simul_repet_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)

#################  Exonic prop according to distance #################
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


################# Sequence conservation human 2 mouse vs distance to promoters #################
obs <- obs[order(obs$ID), ]
simul <- simul[order(simul$ID), ]

obs$align_target <- align_obs[order(align_obs$enh), target_sp]
simul$align_target <- align_simul[order(align_simul$enh), target_sp]

obs_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$align_target)))
obs_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][1])
obs_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][2])

obs_enh_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)))
obs_enh_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][1])
obs_enh_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$align_target)))
simul_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][2])

################# Sequence conservation human 2 macaque vs distance to promoters #################
if(ref_sp == "human"){close_sp="macaque"}else{close_sp="rat"}

obs$align_target <- align_obs[order(align_obs$enh), close_sp]
simul$align_target <- align_simul[order(align_simul$enh), close_sp]

obs_dist_mac <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$align_target)))
obs_dist_mac$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][1])
obs_dist_mac$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][2])

obs_enh_dist_mac <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)))
obs_enh_dist_mac$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][1])
obs_enh_dist_mac$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][2])

simul_dist_mac <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$align_target)))
simul_dist_mac$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][1])
simul_dist_mac$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][2])


####################################### Enhancers datas  ####################################### 
enhancers = c("CAGE", "ENCODE")
if(ref_sp == "human"){enhancers <- c(enhancers, "RoadMap", "GRO_seq")}

list_conserv_enh <- list()
for (enh in enhancers){
  path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, enh, sep="/")
  
  align <- read.table(paste(path_evol,"sequence_conservation", enh, "Alignments_stats_all_species_total_ungapped.txt", sep="/"), header=T)
  obs_stats <- read.table(paste(path_evol,"/sequence_conservation/", enh, "/original_stats.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_evol,"/sequence_conservation/", enh, "/simulated_stats.txt", sep=""), header=T)
  enh_annot <- read.table(paste(path_annot,"/", enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
  
  # Filters
  species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
  align <- align[, c("enh", species)]
  obs_stats$dist_class <- cut(obs_stats$med_dist, breaks=seq(from=25000, to=2025000, by=50000), include.lowest = T)
  simul_stats$dist_class <-cut(simul_stats$med_dist, breaks=seq(from=25000, to=2025000, by=50000), include.lowest = T)
  
  enh_annot_repet <- enh_annot[which(enh_annot$nb_match < 2 & (enh_annot$nb_N/enh_annot$length) < 0.2),]
  
  if (enh == "ENCODE"){
    enh_annot <- enh_annot[which(enh_annot$nb_match < 2 ),]
    align_enhancers_obs <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% obs_stats$enh),]
    align_enhancers_simul <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% simul_stats$enh),]
  }
  
  align_enhancers_dist <- data.frame(matrix(vector(), length(levels(obs_stats$dist_class))))
  
  align_enhancers_dist[[paste0(enh)]] <- sapply(levels(obs_stats$dist_class), function(x) 
    mean(align[which(align$enh %in% enh_annot_repet[,1] & align$enh %in% obs_stats[which(obs_stats$dist_class == x),]$enh),][[target_sp]]))
  
  align_enhancers_dist[[paste0(enh, "_start")]]  <- sapply(levels(obs_stats$dist_class), function(x) 
    t.test(align[which(align$enh %in% enh_annot_repet[,1] & align$enh %in% obs_stats[which(obs_stats$dist_class == x),]$enh),][[target_sp]])[["conf.int"]][1])
  
  align_enhancers_dist[[paste0(enh, "_end")]]  <-sapply(levels(obs_stats$dist_class), function(x) 
    t.test(align[which(align$enh %in% enh_annot_repet[,1] & align$enh %in% obs_stats[which(obs_stats$dist_class == x),]$enh),][[target_sp]])[["conf.int"]][2])
  
  list_conserv_enh <- c(list_conserv_enh, align_enhancers_dist)
  }


# Output
save(align_obs, align_obs_enh, align_simul,
     obs_repet_dist, obs_enh_repet_dist, simul_repet_dist,
     obs_exon_dist, obs_enh_exon_dist, simul_exon_dist,
     obs_dist, obs_enh_dist, simul_dist,
     obs_dist_mac,obs_enh_dist_mac,simul_dist_mac,
     align_enhancers_obs, align_enhancers_simul, list_conserv_enh, 
     file = paste(path, "/Figures/Fig3_", ref_sp, ".Rdata", sep=""))


