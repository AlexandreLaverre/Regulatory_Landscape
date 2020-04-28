ref_sp = "human"
enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")

path <- "/home/laverre/Data/Regulatory_landscape/result/"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_contact <- paste(path, "Supplementary_dataset4_genes_enhancers_contacts", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, "/", sep="/")

################################################# Conserved contact global   #################################################
#pdf(paste(path_evol, "/", ref_sp, "2mouse_contacts_conservation.pdf", sep=""), width=8, height=4.3)
data_global <- c()
conf_low_global <- c()
conf_up_global <- c()
n_total_global <- c()

for (enh in enhancers){
  message(enh)
  all_obs <- read.table(paste(path_contact,"/gene_", enh, "_enhancers_original_interactions.txt", sep=""), header=T, sep="\t")
  all_simul <- read.table(paste(path_contact,"/gene_", enh, "_enhancers_simulated_interactions.txt", sep=""), header=T, sep="\t")
  contact_obs <- read.table(paste(path_evol,"/contact_conservation/", enh, "/human2mouse_original.txt", sep=""), header=T, sep="\t")
  contact_simul <- read.table(paste(path_evol,"/contact_conservation/", enh, "/human2mouse_simulated.txt", sep=""), header=T, sep="\t")
  
  # Select only contacts with no duplicated enhancers
  stats_enh <- read.table(paste(path_annot, enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
  stats_enh <- stats_enh[which(stats_enh$nb_match == 1),]
  all_obs <- all_obs[which(all_obs$enhancer %in% stats_enh$ID),]
  all_simul <- all_simul[which(all_simul$enhancer %in% stats_enh$ID),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% stats_enh$ID),]
  contact_simul <- contact_simul[which(contact_simul$origin_enh %in% stats_enh$ID),]
  
  # Overlap with target enhancer
  if (enh %in% c("CAGE", "ENCODE")){
    overlap_enh <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/", enh, "_lifted_overlap_", enh, "_target.txt", sep=""), header=T, sep="\t")
    enh_tot = nrow(overlap_enh)
    overlap_enh <- overlap_enh[which(overlap_enh$overlap_ID != "NA"),]
    enh_overlap = nrow(overlap_enh) 
    message("Proportion of lifted_enh overlap target enh : ", enh_overlap, " on ", enh_tot, " = ", enh_overlap/enh_tot  )
    
    all_obs <- all_obs[which(all_obs$enhancer %in% overlap_enh$ID),]
    all_simul <- all_simul[which(all_simul$enhancer %in% overlap_enh$ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_enh %in% overlap_enh$ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_enh %in% overlap_enh$ID),]
  }
  
  # Calculate proportion
  mat <- matrix(c(nrow(contact_obs), nrow(contact_simul), nrow(all_obs)-nrow(contact_obs), nrow(all_simul)-nrow(contact_simul)),2)
  
  conf_low_global <- append(conf_low_global, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[1])*100,
                                 (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[1])*100, NA))
  
  conf_up_global <- append(conf_up_global, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[2])*100,
                               (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[2])*100, NA))
  
  data_global <- append(data_global, c( (nrow(contact_obs)/nrow(all_obs))*100, (nrow(contact_simul)/nrow(all_simul))*100, NA))
  
  n_total_global <- append(n_total_global, c(paste0("N = ", nrow(all_obs)), paste0("N = ", nrow(all_simul)), NA))
  
  conserv <- data.frame(data=data_global, conf_low=conf_low_global, conf_up=conf_up_global, n_total=n_total_global)
  assign("conserv_global", conserv)
  
}

################################################# Conserv ~ genomic distance   #################################################
#pdf(paste(path_evol, "/", ref_sp, "2mouse_contacts_conservation_to_distance_difference_obs_simul_1Mb.pdf", sep=""), width=9, height=6)

for (enh in enhancers){
  all_obs <- read.table(paste(path_contact,"/gene_", enh, "_enhancers_original_interactions.txt", sep=""), header=T, sep="\t")
  all_simul <- read.table(paste(path_contact,"/gene_", enh, "_enhancers_simulated_interactions.txt", sep=""), header=T, sep="\t")
  contact_obs <- read.table(paste(path_evol,"/contact_conservation/", enh, "/human2mouse_original.txt", sep=""), header=T, sep="\t")
  contact_simul <- read.table(paste(path_evol,"/contact_conservation/", enh, "/human2mouse_simulated.txt", sep=""), header=T, sep="\t")
  
  # Distance classes
  max_dist = 3000000
  all_obs$class_dist <-cut(all_obs$dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
  all_simul$class_dist <-cut(all_simul$dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
  contact_obs$class_dist <-cut(contact_obs$origin_dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
  contact_simul$class_dist <-cut(contact_simul$origin_dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
  
  class_leg <- c("0",  "500kb",  "1Mb", "1.5Mb", "2Mb", "2.5Mb", "3Mb")
  #class_leg <- c("0",  "1Mb",  "2Mb", "3Mb", "4Mb")
  #class_leg <- c("0",  "250kb",  "500kb", "750kb", "1Mb")
  
  # Select only contacts with no duplicated enhancers
  stats_enh <- read.table(paste(path_annot, enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
  stats_enh <- stats_enh[which(stats_enh$nb_match == 1),]
  all_obs <- all_obs[which(all_obs$enhancer %in% stats_enh$ID),]
  all_simul <- all_simul[which(all_simul$enhancer %in% stats_enh$ID),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% stats_enh$ID),]
  contact_simul <- contact_simul[which(contact_simul$origin_enh %in% stats_enh$ID),]
  
  
  # Calculate proportion
  conserv <- data.frame(result = sapply(levels(all_obs$class_dist), function(x)
    (nrow(contact_obs[which(contact_obs$class_dist == x ),])/nrow(all_obs[which(all_obs$class_dist == x ),]))*100))
  
  conserv$obs_conflow <- sapply(levels(all_obs$class_dist), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_dist == x ),]),
               n=nrow(all_obs[which(all_obs$class_dist == x ),]), p=0.5)$conf.int[1])*100)
  
  conserv$obs_confup <- sapply(levels(all_obs$class_dist), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_dist == x ),]),
               n=nrow(all_obs[which(all_obs$class_dist == x ),]), p=0.5)$conf.int[2])*100)
 
  conserv$simul <- sapply(levels(all_simul$class_dist), function(x)
    (nrow(contact_simul[which(contact_simul$class_dist == x ),])/nrow(all_simul[which(all_simul$class_dist == x ),]))*100)
  
  conserv$simul_conflow <- sapply(levels(all_simul$class_dist), function(x)  
    (prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]),
               n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[1])*100)
  
  conserv$simul_confup <- sapply(levels(all_simul$class_dist), function(x)  
    (prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]),
               n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[2])*100)
  
  assign(paste("conserv_dist", enh, sep = "_"), conserv)

}


################################################# Conserv ~ nb samples #################################################
#pdf(paste(path_evol, "/", ref_sp, "2mouse_contacts_conservation_to_nb_sample.pdf", sep=""), width=8, height=4)

for (enh in enhancers){
  contact_obs <- read.table(paste(path_evol,"/contact_conservation/", enh, "/human2mouse_original.txt", sep=""), header=T, sep="\t")
  contact_simul <- read.table(paste(path_evol,"/contact_conservation/", enh, "/human2mouse_simulated.txt", sep=""), header=T, sep="\t")
  all_obs <- read.table(paste(path_contact,"/gene_", enh, "_enhancers_original_interactions.txt", sep=""), header=T, sep="\t")
  all_simul <- read.table(paste(path_contact,"/gene_", enh, "_enhancers_simulated_interactions.txt", sep=""), header=T, sep="\t")
  
  # Sample classes
  sample_class = c(1, 2, 5, 10, 20)
  all_obs$class_dist <- cut(all_obs$nb_sample, breaks=c(sample_class, max(all_obs$nb_sample)), include.lowest = T)
  all_simul$class_dist <- cut(all_simul$nb_sample, breaks=c(sample_class, max(all_simul$nb_sample)), include.lowest = T)
  contact_obs$class_dist <- cut(contact_obs$nb_sample, breaks=c(sample_class, max(contact_obs$nb_sample)), include.lowest = T)
  contact_simul$class_dist <- cut(contact_simul$nb_sample, breaks=c(sample_class, max(contact_simul$nb_sample)), include.lowest = T)
  
  # Select only contacts with no duplicated enhancers
  stats_enh <- read.table(paste(path_annot, enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
  stats_enh <- stats_enh[which(stats_enh$nb_match == 1),]
  all_obs <- all_obs[which(all_obs$enhancer %in% stats_enh$ID),]
  all_simul <- all_simul[which(all_simul$enhancer %in% stats_enh$ID),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% stats_enh$ID),]
  contact_simul <- contact_simul[which(contact_simul$origin_enh %in% stats_enh$ID),]
  
  # Overlap with target enhancer
  if (enh %in% c("CAGE", "ENCODE")){
    overlap_enh <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/", enh, "_lifted_overlap_", enh, "_target.txt", sep=""), header=T, sep="\t")
    enh_tot = nrow(overlap_enh)
    overlap_enh <- overlap_enh[which(overlap_enh$overlap_ID != "NA"),]
    enh_overlap = nrow(overlap_enh) 
    message("Proportion of lifted_enh overlap target enh : ", enh_overlap, " on ", enh_tot, " = ", enh_overlap/enh_tot  )
    
    all_obs <- all_obs[which(all_obs$enhancer %in% overlap_enh$ID),]
    all_simul <- all_simul[which(all_simul$enhancer %in% overlap_enh$ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_enh %in% overlap_enh$ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_enh %in% overlap_enh$ID),]
  }
  
  # Calculate proportion
  conserv_obs <- data.frame(result = sapply(levels(all_obs$class_dist), function(x)
    (nrow(contact_obs[which(contact_obs$class_dist == x),])/(nrow(all_obs[which(all_obs$class_dist == x ),])))*100))
  
  conserv_obs$conflow <- sapply(levels(all_obs$class_dist), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_dist == x ),]),
               n=nrow(all_obs[which(all_obs$class_dist == x ),]), p=0.5)$conf.int[1])*100)
  
  conserv_obs$confup <- sapply(levels(all_obs$class_dist), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_dist == x ),]),
               n=nrow(all_obs[which(all_obs$class_dist == x ),]), p=0.5)$conf.int[2])*100)
  
  conserv_obs$data <- "obs"
  conserv_obs$class <-  c("1", "2-5", "6-10", "11-20", "21-26")
  conserv_obs$count <- sapply(levels(all_obs$class_dist), function(x) paste("n=", nrow(all_obs[which(all_obs$class_dist == x ),]), sep=""))
  
  conserv_simul <- data.frame(result = sapply(levels(all_simul$class_dist), function(x)
    (nrow(contact_simul[which(contact_simul$class_dist == x ),])/(nrow(all_simul[which(all_simul$class_dist == x ),])))*100))
  
  conserv_simul$conflow <- sapply(levels(all_simul$class_dist), function(x) 
    tryCatch((prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]),
               n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[1])*100,error=function(e) 0))

  conserv_simul$confup <- sapply(levels(all_simul$class_dist), function(x)  
    tryCatch((prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]),
                        n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[2])*100,error=function(e) 0))
  
  conserv_simul$data <- "simul"
  conserv_simul$class <- c("1", "2-5", "6-10", "11-20", "21-26")
  conserv_simul$count <- sapply(levels(all_simul$class_dist), function(x) paste("n=", nrow(all_simul[which(all_simul$class_dist == x ),]), sep=""))
  conserv <- rbind(conserv_obs, conserv_simul)
  conserv$class <- ordered(conserv$class, levels = c("1", "2-5", "6-10", "11-20", "21-26"))
  
  assign(paste("conserv_sample", enh, sep = "_"), conserv)
  
}

################################################# Conserv in similar samples #################################################
for (enh in enhancers){
  message(enh)
  data <- c()
  conf_low <- c()
  conf_up <- c()
  p_test <- c()
  n_total <- c()
  
  all_obs <- read.table(paste(path_contact,"/gene_", enh, "_enhancers_original_interactions.txt", sep=""), header=T, sep="\t")
  all_simul <- read.table(paste(path_contact,"/gene_", enh, "_enhancers_simulated_interactions.txt", sep=""), header=T, sep="\t")
  contact_obs <- read.table(paste(path_evol,"/contact_conservation/", enh, "/human2mouse_original.txt", sep=""), header=T, sep="\t")
  contact_simul <- read.table(paste(path_evol,"/contact_conservation/", enh, "/human2mouse_simulated.txt", sep=""), header=T, sep="\t")
  
  # Select only contacts with no duplicated enhancers
  stats_enh <- read.table(paste(path_annot, enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
  stats_enh <- stats_enh[which(stats_enh$nb_match == 1),]
  all_obs <- all_obs[which(all_obs$enhancer %in% stats_enh$ID),]
  all_simul <- all_simul[which(all_simul$enhancer %in% stats_enh$ID),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% stats_enh$ID),]
  contact_simul <- contact_simul[which(contact_simul$origin_enh %in% stats_enh$ID),]
  
  # Overlap with target enhancer
  if (enh %in% c("CAGE", "ENCODE")){
    overlap_enh <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/", enh, "_lifted_overlap_", enh, "_target.txt", sep=""), header=T, sep="\t")
    enh_tot = nrow(overlap_enh)
    overlap_enh <- overlap_enh[which(overlap_enh$overlap_ID != "NA"),]
    enh_overlap = nrow(overlap_enh) 
    message("Proportion of lifted_enh overlap target enh : ", enh_overlap, " on ", enh_tot, " = ", enh_overlap/enh_tot  )
    
    all_obs <- all_obs[which(all_obs$enhancer %in% overlap_enh$ID),]
    all_simul <- all_simul[which(all_simul$enhancer %in% overlap_enh$ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_enh %in% overlap_enh$ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_enh %in% overlap_enh$ID),]
  }
  
  
  ### Select relative cell type
  # Pre-adipocytes
  all_obs_adip <- all_obs[which(!is.na(all_obs$pre_adipo)),]
  all_simul_adip <- all_simul[which(!is.na(all_simul$pre_adipo)),]
  contact_obs_adip <- contact_obs[which(contact_obs$pre_adipo > 0 & (contact_obs$preadip_D0 > 0 | contact_obs$preadip_D2 > 0 | contact_obs$preadip_4H > 0 )),]
  contact_simul_adip <- contact_simul[which(contact_simul$pre_adipo > 0 & (contact_simul$preadip_D0 > 0 | contact_simul$preadip_D2 > 0 | contact_simul$preadip_4H > 0 )),]
  
  # ESC
  all_obs_ESC <- all_obs[which(!is.na(all_obs$hESC)),]
  all_simul_ESC <- all_simul[which(!is.na(all_simul$hESC)),]
  contact_obs_ESC <- contact_obs[which(contact_obs$hESC > 0 & (contact_obs$ESC > 0 | contact_obs$ESC_18 > 0 | contact_obs$ESC_NKO > 0 | contact_obs$ESC_wild > 0 )),]
  contact_simul_ESC <- contact_simul[which(contact_simul$hESC > 0 & (contact_simul$ESC > 0 | contact_simul$ESC_18 > 0 | contact_simul$ESC_NKO > 0 | contact_simul$ESC_wild > 0)),]
  
  # Bcell
  all_obs_Bcell <- all_obs[which(!is.na(all_obs$Bcell) | !is.na(all_obs$TB)),]
  all_simul_Bcell <- all_simul[which(!is.na(all_simul$Bcell) | !is.na(all_simul$TB)),]
  contact_obs_Bcell <- contact_obs[which((contact_obs$Bcell > 0 | contact_obs$TB > 0) & (contact_obs$preB_aged > 0 | contact_obs$preB_young > 0)),]
  contact_simul_Bcell <- contact_simul[which((contact_simul$Bcell > 0 | contact_simul$TB > 0) & (contact_simul$preB_aged > 0 | contact_simul$preB_young > 0)),]
  
  for (cell in c("adip", "ESC", "Bcell")){
    # Calculate proportion
    all_obs = get(paste("all_obs", cell, sep="_"))
    all_simul = get(paste("all_simul", cell, sep="_"))
    contact_obs = get(paste("contact_obs", cell, sep="_"))
    contact_simul = get(paste("contact_simul", cell, sep="_"))
    
    mat <- matrix(c(nrow(contact_obs), nrow(contact_simul), nrow(all_obs)-nrow(contact_obs), nrow(all_simul)-nrow(contact_simul)),2)
    
    conf_low <- append(conf_low, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[1])*100,
                                   (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[1])*100, NA))
    
    conf_up <- append(conf_up, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[2])*100,
                                 (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[2])*100, NA))
    
    data <- append(data, c( (nrow(contact_obs)/nrow(all_obs))*100, (nrow(contact_simul)/nrow(all_simul))*100, NA))
    
    n_total <- append(n_total, c(paste0("N = ", nrow(all_obs)), paste0("N = ", nrow(all_simul)), NA))
    
  }

  conserv <- data.frame(data=data, conf_low=conf_low, conf_up=conf_up, n_total=n_total)
  assign(paste("conserv_similar_sample", enh, sep = "_"), conserv)
  
}

save(conserv_global,
     conserv_dist_CAGE, conserv_dist_ENCODE, conserv_dist_RoadMap, conserv_dist_GRO_seq,
     conserv_sample_CAGE, conserv_sample_ENCODE, conserv_sample_RoadMap, conserv_sample_GRO_seq,
     conserv_similar_sample_CAGE, conserv_similar_sample_ENCODE, conserv_similar_sample_RoadMap, conserv_similar_sample_GRO_seq,
     file = paste(path, "Main_figures/Fig5_human_corrected.Rdata", sep=""))
