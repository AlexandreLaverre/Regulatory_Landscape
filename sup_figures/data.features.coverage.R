#########################################################################################################################

library(data.table)

options(stringsAsFactors = FALSE)

source("../main_figures/parameters.R")

library(bootBCa, lib=pathRlibs)

set.seed(19)

load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))

#########################################################################################################################

for(ref_sp in c("human", "mouse")){
  feat_prop_dist <- list()
  
  enhancers = enhancer.datasets[[ref_sp]]
  features = c("exon", "repeat", "GC", "genes")

  ## read fragment statistics, already filtered
  obs=fragment.statistics[[ref_sp]][["original"]]
  simul=fragment.statistics[[ref_sp]][["simulated"]]
  
  ## compute percentage of length covered by other features
  obs$exon_pclen=obs$exon_bp*100/obs$length
  simul$exon_pclen=simul$exon_bp*100/simul$length
  
  obs$repeat_pclen=obs$repeat_bp*100/obs$length
  simul$repeat_pclen=simul$repeat_bp*100/simul$length
  
  obs$GC_pclen=100*obs$GC_content
  simul$GC_pclen=100*simul$GC_content
  
  obs$genes_pclen=obs$nb_genes_500kb
  simul$genes_pclen=simul$nb_genes_500kb

  ############################### features coverage according to distance from promoters #################
  
  obs$dist_class <-cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  simul$dist_class <- cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  
  obs_feat_dist <- list()
  simul_feat_dist <- list()
  
  for (feat in features){
    print(paste("confidence intervals fragments", feat))
    print("observed")
    BC.obs=tapply(obs[, paste0(feat, "_pclen")], obs$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    obs_feat_dist[[feat]] <- unlist(lapply(BC.obs, function(x) x[3]))
    obs_feat_dist[[paste0(feat, "_conflow")]] <-  unlist(lapply(BC.obs, function(x) x[4]))
    obs_feat_dist[[paste0(feat, "_confup")]] <-  unlist(lapply(BC.obs, function(x) x[5]))

    print("simulated")
    BC.simul=tapply(simul[, paste0(feat, "_pclen")], simul$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    simul_feat_dist[[feat]] <- unlist(lapply(BC.simul, function(x) x[3]))
    simul_feat_dist[[paste0(feat, "_conflow")]] <-  unlist(lapply(BC.simul, function(x) x[4]))
    simul_feat_dist[[paste0(feat, "_confup")]] <-  unlist(lapply(BC.simul, function(x) x[5]))
  }
  
  feat_prop_dist[["fragment"]] <- list(obs=obs_feat_dist, simul=simul_feat_dist)
  
  #########################################################################################################################
  #########################################################################################################################
############################### features coverage in enhancers according to distance from promoters #################
  
  for (enh in enhancers){

    obs=enhancer.statistics[[ref_sp]][[enh]][["original"]]
    simul=enhancer.statistics[[ref_sp]][[enh]][["simulated"]]
   
    ## compute percentage of length covered by other features
    obs$exon_pclen=obs$exon_bp*100/obs$length
    simul$exon_pclen=simul$exon_bp*100/simul$length
    
    obs$repeat_pclen=obs$repeat_bp*100/obs$length
    simul$repeat_pclen=simul$repeat_bp*100/simul$length
    
    obs$GC_pclen=100*obs$GC_content
    simul$GC_pclen=100*simul$GC_content

    obs$genes_pclen=obs$nb_genes_500kb
    simul$genes_pclen=simul$nb_genes_500kb
        
    ############################### features coverage according to distance from promoters #################
    obs$dist_class <-cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    simul$dist_class <- cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    
    obs_feat_dist <- list()
    simul_feat_dist <- list()
    
    for (feat in features){

      print(paste("confidence intervals",enh, feat))
      print("observed")
      BC.obs=tapply(obs[, paste0(feat, "_pclen")], obs$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
      obs_feat_dist[[feat]] <- unlist(lapply(BC.obs, function(x) x[3]))
      obs_feat_dist[[paste0(feat, "_conflow")]] <-  unlist(lapply(BC.obs, function(x) x[4]))
      obs_feat_dist[[paste0(feat, "_confup")]] <-  unlist(lapply(BC.obs, function(x) x[5]))

      print("simulated")
      BC.simul=tapply(simul[, paste0(feat, "_pclen")], simul$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
      simul_feat_dist[[feat]] <- unlist(lapply(BC.simul, function(x) x[3]))
      simul_feat_dist[[paste0(feat, "_conflow")]] <-  unlist(lapply(BC.simul, function(x) x[4]))
      simul_feat_dist[[paste0(feat, "_confup")]] <-  unlist(lapply(BC.simul, function(x) x[5]))
    }
    
    feat_prop_dist[[enh]] <- list(obs=obs_feat_dist, simul=simul_feat_dist)
    
  }
  
  ################################################# Save RData ################################################# 
  
  save(feat_prop_dist, file = paste(pathFigures, "RData/data.features.coverage.", ref_sp, ".RData", sep=""))
  
  ############################################################################################################## 
  
}

############################################################################################################## 
