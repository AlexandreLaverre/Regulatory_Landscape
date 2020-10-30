################################################################################################################################################

library(data.table)

options(stringsAsFactors = FALSE)

source("parameters.R") ## paths are defined based on the user name

################################################################################################################################################

for(ref_sp in c("human", "mouse")){
  print(ref_sp)
  
  if(ref_sp == "human"){
    target_sp = "mouse"
    closest_sp="macaque"
  }else{
    target_sp = "human"
    closest_sp="rat"
  }

  enhancers=enhancer.datasets[[ref_sp]]
  
  species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
  
  path_evol <- paste(pathFinalData, "SupplementaryDataset7/", ref_sp, "/", sep="")
  path_annot <- paste(pathFinalData, "SupplementaryDataset4/", ref_sp, "/", sep="")
  
  ####  statistics for the contacted restriction fragments
  obs <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp, "/statistics_contacted_sequence_original.txt", sep=""), header=T)
  simul <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp,"/statistics_contacted_sequence_simulated.txt", sep=""), header=T)
  
  obs$ID <-  do.call(paste,c(obs[c("chr","start","end")],sep=":"))
  simul$ID <-  do.call(paste,c(simul[c("chr","start","end")],sep=":"))

  rownames(obs) <- obs$ID
  rownames(simul) <- simul$ID
  
  # unbaited only
  obs <- obs[which(obs$baited == "unbaited"),]
  simul <- simul[which(simul$baited == "unbaited"),]

  ## remove duplicated sequences
  
  obs <- obs[which(obs$BLAT_match < 2),]
  simul <- simul[which(simul$BLAT_match < 2),]

  #### alignment score vs all species, for restriction fragments 
  frag_align <- read.table(paste(path_evol,"sequence_conservation/restriction_fragments/Alignments_stats_all_species_nonexonic_ungapped.txt", sep="/"), header=T)
    
  ID=paste("ID", ref_sp, sep=".")
  
  frag_align_obs <- frag_align[which(frag_align[,ID] %in% obs$ID), c(ID, species)]
  frag_align_simul <- frag_align[which(frag_align[,ID] %in% simul$ID), c(ID, species) ]

  frag_align_obs$median_dist <- obs[frag_align_obs[,ID],"median_dist"]
  frag_align_simul$median_dist <- simul[frag_align_simul[,ID],"median_dist"]

  frag_align_obs$dist_class <- cut(frag_align_obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  frag_align_simul$dist_class <-cut(frag_align_simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      
  ## alignment for enhancers
  
  list_align_enh <- list()
  
  for (enh in enhancers){
    print(enh)
    
    enh_align <- read.table(paste(path_evol,"sequence_conservation/enhancers/", enh, "Alignments_stats_all_species_nonexonic_ungapped.txt", sep="/"), header=T)
    enh_obs_stats <- read.table(paste(path_annot, enh, "/statistics_contacted_enhancers_original.txt", sep=""), header=T)
    enh_simul_stats <- read.table(paste(path_annot, enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), header=T)
    
    enh_obs_stats$enh <-  do.call(paste,c(enh_obs_stats[c("chr","start","end")],sep=":"))
    enh_simul_stats$enh <-  do.call(paste,c(enh_simul_stats[c("chr","start","end")],sep=":"))

    rownames(enh_obs_stats) <- enh_obs_stats[,"enh"]
    rownames(enh_simul_stats) <- enh_simul_stats[,"enh"]

    ## distance class
    enh_obs_stats$dist_class <- cut(enh_obs_stats$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    enh_simul_stats$dist_class <-cut(enh_simul_stats$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    
    
    ## Select unduplicated and with repeat_part < 1%
    enh_obs_stats <- enh_obs_stats[which(enh_obs_stats$BLAT_match < 2 & (enh_obs_stats$repeat_bp/enh_obs_stats$length) < 0.01),] 
    enh_simul_stats <- enh_simul_stats[which(enh_simul_stats$BLAT_match < 2 & (enh_simul_stats$repeat_bp/enh_simul_stats$length) < 0.01),]
    
    ## select species
    enh_align <- enh_align[, c("enh", species)]
    
    enh_align_obs <- enh_align[which(enh_align$enh %in% enh_obs_stats$enh),]
    enh_align_simul<- enh_align[which(enh_align$enh %in% enh_simul_stats$enh),]

    enh_align_obs$median_dist <- enh_obs_stats[enh_align_obs$enh, "median_dist"]
    enh_align_obs$dist_class <- enh_obs_stats[enh_align_obs$enh, "dist_class"]

    enh_align_simul$median_dist <- enh_simul_stats[enh_align_simul$enh, "median_dist"]
    enh_align_simul$dist_class <- enh_simul_stats[enh_align_simul$enh, "dist_class"]

    ## store results
    list_align_enh[[enh]] <- list("enh_align_obs", "enh_align_simul")
  }

  ## save results

  save(list=c("frag_align_obs", "frag_align_simul", "list_align_enh"), file=paste(pathFigures, "RData/data.sequence.conservation.",ref_sp,".Rdata", sep=""))
  
}

################################################################################################################################################
