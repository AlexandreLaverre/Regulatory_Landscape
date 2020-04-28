ref_sp = "human"
target_sp = "mouse" 
enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")
col <- c("red", "navy", "forestgreen", "orange")

path <- "/home/laverre/Data/Regulatory_landscape/result"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, sep="/")

enh <- "ENCODE"


# BLAT
enh_annot <- read.table(paste(path_annot,"/", enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
enh_annot <- enh_annot[which(enh_annot$nb_match == 1), ]
# conserv
align <- read.table(paste(path_evol,"enhancers_conservation", enh, "Alignments_stats_all_species_identity.txt", sep="/"), header=T)
align <- align[which(align$mouse > 0.1 & align$enh %in% enh_annot$ID),]
# overlap
overlap_enh <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/", enh, "_lifted_overlap_", enh, "_target.txt", sep=""), header=T, sep="\t")
overlap_enh <- overlap_enh[which(overlap_enh$overlap_ID != "NA" & overlap_enh$ID %in% align$enh),]
# enh in contact
obs_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/original_stats.txt", sep=""), header=T)
simul_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/simulated_stats.txt", sep=""), header=T)


# BLAT
tot_obs <- nrow(obs_stats)
tot_obs
obs_stats <- obs_stats[which(obs_stats$enh %in% enh_annot$ID),]
blat_obs <- nrow(obs_stats)
blat_obs/tot_obs

tot_simul <- nrow(simul_stats)
tot_simul
simul_stats <- simul_stats[which(simul_stats$enh %in% enh_annot$ID),]
blat_simul <- nrow(simul_stats)
blat_simul/tot_simul

# conserv
obs_stats <- obs_stats[which(obs_stats$enh %in% align$enh),]
conserv_obs <- nrow(obs_stats)
conserv_obs/blat_obs

simul_stats <- simul_stats[which(simul_stats$enh %in% align$enh),]
conserv_simul <- nrow(simul_stats)
conserv_simul/blat_simul

# overlap
obs_stats <- obs_stats[which(obs_stats$enh %in% overlap_enh$ID),]
overlap_obs <- nrow(obs_stats)
overlap_obs/conserv_obs

simul_stats <- simul_stats[which(simul_stats$enh %in% overlap_enh$ID),]
overlap_simul <- nrow(simul_stats)
overlap_simul/conserv_simul
