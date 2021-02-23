#########################################################################################################################

source("parameters.R")

library(data.table)
library(bootBCa, lib=pathRlibs)

set.seed(19)

options(stringsAsFactors = FALSE)

load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))

#########################################################################################################################

for(ref_sp in c("human", "mouse")){
  
  enhancers = enhancer.datasets[[ref_sp]]
  
  obs <- fragment.statistics[[ref_sp]][["original"]]
  simul <- fragment.statistics[[ref_sp]][["simulated"]]

  ## compute percentage of length covered by enhancers
  
  for(enh in enhancers){
    obs[,paste0(enh,"_pclen")]=as.numeric(obs[,paste0(enh, "_bp")])*100/obs$length
    simul[,paste0(enh,"_pclen")]=as.numeric(simul[,paste0(enh, "_bp")])*100/simul$length
  }
  
  overlap = c(length(which(as.numeric(obs[,"ENCODE_bp"]) > 0)), length(which(as.numeric(simul[,"ENCODE_bp"]) > 0)))
  non.overlap = c(length(which(as.numeric(obs[,"ENCODE_bp"]) == 0)), length(which(as.numeric(simul[,"ENCODE_bp"]) == 0)))
  
  print(paste(ref_sp, "observed restriction fragment overlap with at least one ENCODE enhancer", round((overlap[1]*100/nrow(obs)),2), "%"))
  print(paste(ref_sp, "simulated restriction fragment overlap with at least one ENCODE enhancer", round((overlap[2]*100/nrow(simul)),2), "%"))
  print(chisq.test(matrix(c(overlap, non.overlap), nrow=2)))
  
  ###################### Proportions of contacted sequences covered by enhancers ########################
  
  data <- c()
  conf_up <- c()
  conf_low <- c()
  id <- c()
    
  ### Proportion of the sequences covered by enhancers

  for (enh in enhancers){
    print(paste("confidence intervals, total",enh))
    
    x <- BCa(obs[,paste0(enh, "_pclen")], delta=NA, M=100, theta=mean, na.rm=T)
    data <- c(data, x[3])
    conf_up <- c(conf_up, x[5])
    conf_low <- c(conf_low, x[4])
    id <- c(id, paste0(enh,":obs"))
    
    x <- BCa(obs[,paste0(enh, "_pclen")], delta=NA, M=100, theta=mean, na.rm=T)
    data <- c(data, x[3])
    conf_up <- c(conf_up, x[5])
    conf_low <- c(conf_low, x[4])
    id <- c(id, paste0(enh,":sim"))
  }
  
  enh_prop <- data.frame(id=id, data=data, conf_up=conf_up, conf_low=conf_low)
  
  ############################### enhancer coverage according to distance from promoters #################
  
  obs$dist_class <-cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  simul$dist_class <- cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  
  obs_enh_dist <- list()
  simul_enh_dist <- list()
  
  for (enh in enhancers){
    print(paste("confidence intervals by distance",enh))
    
    BC.obs=tapply(obs[, paste0(enh, "_pclen")], obs$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    obs_enh_dist[[enh]] <- unlist(lapply(BC.obs, function(x) x[3]))
    obs_enh_dist[[paste0(enh, "_conflow")]] <- unlist(lapply(BC.obs, function(x) x[4]))
    obs_enh_dist[[paste0(enh, "_confup")]] <- unlist(lapply(BC.obs, function(x) x[5]))

    BC.sim=tapply(simul[, paste0(enh, "_pclen")], simul$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    simul_enh_dist[[enh]] <- unlist(lapply(BC.sim, function(x) x[3]))
    simul_enh_dist[[paste0(enh, "_conflow")]] <- unlist(lapply(BC.sim, function(x) x[4]))
    simul_enh_dist[[paste0(enh, "_confup")]] <- unlist(lapply(BC.sim, function(x) x[5]))
  }
  
  enh_prop_dist <- list(obs=obs_enh_dist, simul=simul_enh_dist)

  ###############################  proportion of sequence covered by enhancers by cell type  ############
  
  obs$nb_cell <- as.factor(obs$nb_cell)
  simul$nb_cell <- as.factor(simul$nb_cell)
  
  obs_enh_cell <- list()
  simul_enh_cell <- list()
  
  for (enh in enhancers){
    print(paste("confidence intervals by cell type",enh))
    
    BC.obs=tapply(obs[, paste0(enh, "_pclen")], obs$nb_cell, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    obs_enh_cell[[enh]] <- unlist(lapply(BC.obs, function(x) x[3]))
    obs_enh_cell[[paste0(enh, "_conflow")]] <- unlist(lapply(BC.obs, function(x) x[4]))
    obs_enh_cell[[paste0(enh, "_confup")]] <- unlist(lapply(BC.obs, function(x) x[5]))

    BC.sim=tapply(simul[, paste0(enh, "_pclen")], simul$nb_cell, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    simul_enh_cell[[enh]] <- unlist(lapply(BC.sim, function(x) x[3]))
    simul_enh_cell[[paste0(enh, "_conflow")]] <- unlist(lapply(BC.sim, function(x) x[4]))
    simul_enh_cell[[paste0(enh, "_confup")]] <- unlist(lapply(BC.sim, function(x) x[5]))
  }
  
  enh_prop_nb_cell <- list(obs=obs_enh_cell, simul=simul_enh_cell)
  
  ################################################# Save RData ################################################# 
  
  save(enh_prop, enh_prop_nb_cell, enh_prop_dist, file = paste(pathFigures, "RData/data.enhancer.coverage.", ref_sp, ".RData", sep=""))

  ############################################################################################################## 

}

############################################################################################################## 
