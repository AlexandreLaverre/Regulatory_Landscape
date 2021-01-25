#########################################################################################################################

library(data.table)

options(stringsAsFactors = FALSE)

source("../main_figures/parameters.R") 

#########################################################################################################################

for(ref_sp in c("human", "mouse")){
  feat_prop_dist <- list()
  
  enhancers = enhancer.datasets[[ref_sp]]
  features = c("all_exon", "repeat", "GC")
  
  obs <- fread(paste(pathFinalData, "SupplementaryDataset5/", ref_sp, "/statistics_contacted_sequence_original.txt", sep=""), header=T)
  simul <- fread(paste(pathFinalData, "SupplementaryDataset5/", ref_sp,"/statistics_contacted_sequence_simulated.txt", sep=""), header=T)
  
  class(obs) <- "data.frame"
  class(simul) <- "data.frame"
  
  obs <- obs[which(obs$baited == "unbaited" & obs$BLAT_match == 1 ),]
  simul <- simul[which(simul$baited == "unbaited" & simul$BLAT_match == 1),]
  
  ## compute percentage of length covered by other features
  obs$all_exon_pclen=obs$all_exon_bp*100/obs$length
  simul$all_exon_pclen=simul$all_exon_bp*100/simul$length
  
  obs$repeat_pclen=obs$repet_noexon_bp*100/obs$length
  simul$repeat_pclen=simul$repet_noexon_bp*100/simul$length
  
  obs$GC_pclen=obs$GC_bp*100/(obs$length-obs$repeat_bp)
  simul$GC_pclen=simul$GC_bp*100/(simul$length-simul$repeat_bp)

  ############################### features coverage according to distance from promoters #################
  
  obs$dist_class <-cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  simul$dist_class <- cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  
  obs_feat_dist <- list()
  simul_feat_dist <- list()
  
  for (feat in features){
    obs_feat_dist[[feat]] <- tapply(obs[, paste0(feat, "_pclen")], obs$dist_class, mean, na.rm=T)
    obs_feat_dist[[paste0(feat, "_conflow")]] <- tapply(obs[, paste0(feat, "_pclen")], obs$dist_class, function(x) t.test(x)[["conf.int"]][1])
    obs_feat_dist[[paste0(feat, "_confup")]] <- tapply(obs[, paste0(feat, "_pclen")], obs$dist_class, function(x) t.test(x)[["conf.int"]][2])
    
    simul_feat_dist[[feat]] <- tapply(simul[, paste0(feat, "_pclen")], simul$dist_class, mean, na.rm=T)
    simul_feat_dist[[paste0(feat, "_conflow")]] <- tapply(simul[, paste0(feat, "_pclen")], simul$dist_class, function(x) t.test(x)[["conf.int"]][1])
    simul_feat_dist[[paste0(feat, "_confup")]] <- tapply(simul[, paste0(feat, "_pclen")], simul$dist_class, function(x) t.test(x)[["conf.int"]][2])
  }
  
  feat_prop_dist[["fragment"]] <- list(obs=obs_feat_dist, simul=simul_feat_dist)
  
  #########################################################################################################################
  #########################################################################################################################
  ############################### features coverage in enhancers according to distance from promoters #################
  for (enh in enhancers){
    
    obs <- fread(paste(pathFinalData, "SupplementaryDataset4/", ref_sp, "/", enh, "/statistics_contacted_enhancers_original.txt", sep=""), header=T)
    simul <- fread(paste(pathFinalData, "SupplementaryDataset4/", ref_sp,"/", enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), header=T)
    
    obs <- obs[which(obs$BLAT_match == 1),]
    simul <- simul[which(simul$BLAT_match == 1),]
    
    class(obs) <- "data.frame"
    class(simul) <- "data.frame"
    
    ## compute percentage of length covered by other features
    obs$all_exon_pclen=obs$all_exon_bp*100/obs$length
    simul$all_exon_pclen=simul$all_exon_bp*100/simul$length
    
    obs$repeat_pclen=obs$repeat_bp*100/obs$length
    simul$repeat_pclen=simul$repeat_bp*100/simul$length
    
    obs$GC_pclen=obs$GC_bp*100/(obs$length-obs$repeat_bp)
    simul$GC_pclen=simul$GC_bp*100/(simul$length-simul$repeat_bp)
    
    ############################### features coverage according to distance from promoters #################
    obs$dist_class <-cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    simul$dist_class <- cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    
    obs_feat_dist <- list()
    simul_feat_dist <- list()
    
    for (feat in features){
      obs_feat_dist[[feat]] <- tapply(obs[, paste0(feat, "_pclen")], obs$dist_class, mean, na.rm=T)
      obs_feat_dist[[paste0(feat, "_conflow")]] <- tapply(obs[, paste0(feat, "_pclen")], obs$dist_class, function(x) 
        tryCatch(t.test(x)[["conf.int"]][1], error = function(e) 0))
      obs_feat_dist[[paste0(feat, "_confup")]] <- tapply(obs[, paste0(feat, "_pclen")], obs$dist_class, function(x)
        tryCatch(t.test(x)[["conf.int"]][2], error = function(e) 0))
      
      simul_feat_dist[[feat]] <- tapply(simul[, paste0(feat, "_pclen")], simul$dist_class, mean, na.rm=T)
      simul_feat_dist[[paste0(feat, "_conflow")]] <- tapply(simul[, paste0(feat, "_pclen")], simul$dist_class, function(x)
        tryCatch(t.test(x)[["conf.int"]][1], error = function(e) 0))
      simul_feat_dist[[paste0(feat, "_confup")]] <- tapply(simul[, paste0(feat, "_pclen")], simul$dist_class, function(x)
        tryCatch(t.test(x)[["conf.int"]][2], error = function(e) 0))
    }
    
    feat_prop_dist[[enh]] <- list(obs=obs_feat_dist, simul=simul_feat_dist)
    
  }
  
  ################################################# Save RData ################################################# 
  
  save(feat_prop_dist, file = paste(pathFigures, "RData/data.features.coverage.", ref_sp, ".Rdata", sep=""))
  
  ############################################################################################################## 
  
}

############################################################################################################## 
