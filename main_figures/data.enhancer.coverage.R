#########################################################################################################################

library(data.table)

options(stringsAsFactors = FALSE)

source("parameters.R") 

#########################################################################################################################

for(ref_sp in c("human", "mouse")){
  
  enhancers = enhancer.datasets[[ref_sp]]
  
  obs <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp, "/statistics_contacted_sequence_original.txt", sep=""), header=T)
  simul <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp,"/statistics_contacted_sequence_simulated.txt", sep=""), header=T)
  
  ## obs <- obs[which(obs$baited == "unbaited" & obs$BLAT_match == 1),]
  ## simul <- simul[which(simul$baited == "unbaited" & simul$BLAT_match == 1),]

  obs <- obs[which(obs$baited == "unbaited"),]
  simul <- simul[which(simul$baited == "unbaited"),]


  ## compute percentage of length covered by enhancers
  
  for(enh in enhancers){
    obs[,paste0(enh,"_pclen")]=obs[,paste0(enh, "_bp")]*100/obs$length
    simul[,paste0(enh,"_pclen")]=simul[,paste0(enh, "_bp")]*100/obs$length
  }
  
  ###################### Proportions of contacted sequences covered by enhancers ########################
  
  data <- c()
  conf_up <- c()
  conf_low <- c()
  id <- c()
    
  ### Proportion of the sequences covered by enhancers

  for (enh in enhancers){
    x <- t.test(obs[,paste0(enh, "_pclen")])
    data <- c(data, x$estimate)
    conf_up <- c(conf_up, x$conf.int[1])
    conf_low <- c(conf_low, x$conf.int[2])
    id <- c(id, paste0(enh,":obs"))
    
    x <- t.test(simul[,paste0(enh, "_pclen")])
    data <- c(data, x$estimate)
    conf_up <- c(conf_up, x$conf.int[1])
    conf_low <- c(conf_low, x$conf.int[2])
    id <- c(id, paste0(enh,":sim"))
  }
  
  enh_prop <- data.frame(id=id, data=data, conf_up=conf_up, conf_low=conf_low)
  
  ############################### enhancer coverage according to distance from promoters #################
  
  obs$dist_class <-cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  simul$dist_class <- cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
  
  obs_enh_dist <- list()
  simul_enh_dist <- list()
  
  for (enh in enhancers){
    obs_enh_dist[[enh]] <- tapply(obs[, paste0(enh, "_pclen")], obs$dist_class, mean)
    obs_enh_dist[[paste0(enh, "_conflow")]] <- tapply(obs[, paste0(enh, "_pclen")], obs$dist_class, function(x) t.test(x)[["conf.int"]][1])
    obs_enh_dist[[paste0(enh, "_confup")]] <- tapply(obs[, paste0(enh, "_pclen")], obs$dist_class, function(x) t.test(x)[["conf.int"]][2])

    simul_enh_dist[[enh]] <- tapply(simul[, paste0(enh, "_pclen")], simul$dist_class, mean)
    simul_enh_dist[[paste0(enh, "_conflow")]] <- tapply(simul[, paste0(enh, "_pclen")], simul$dist_class, function(x) t.test(x)[["conf.int"]][1])
    simul_enh_dist[[paste0(enh, "_confup")]] <- tapply(simul[, paste0(enh, "_pclen")], simul$dist_class, function(x) t.test(x)[["conf.int"]][2])
  }
  
  enh_prop_dist <- list(obs=obs_enh_dist, simul=simul_enh_dist)

  ###############################  proportion of sequence covered by enhancers by cell type  ############
  
  obs$nb_cell <- as.factor(obs$nb_cell)
  simul$nb_cell <- as.factor(simul$nb_cell)
  
  obs_enh_cell <- list()
  simul_enh_cell <- list()
  
  for (enh in enhancers){
    obs_enh_cell[[enh]] <- tapply(obs[, paste0(enh, "_pclen")], obs$nb_cell, mean)
    obs_enh_cell[[paste0(enh, "_conflow")]] <- tapply(obs[, paste0(enh, "_pclen")], obs$nb_cell, function(x) tryCatch(t.test(x)[["conf.int"]][1], error=function(e) 0))
    obs_enh_cell[[paste0(enh, "_confup")]] <- tapply(obs[, paste0(enh, "_pclen")], obs$nb_cell, function(x) tryCatch(t.test(x)[["conf.int"]][2], error=function(e) 0))


    simul_enh_cell[[enh]] <- tapply(simul[, paste0(enh, "_pclen")], simul$nb_cell, mean)
    simul_enh_cell[[paste0(enh, "_conflow")]] <- tapply(simul[, paste0(enh, "_pclen")], simul$nb_cell, function(x) tryCatch(t.test(x)[["conf.int"]][1], error=function(e) 0))
    simul_enh_cell[[paste0(enh, "_confup")]] <- tapply(simul[, paste0(enh, "_pclen")], simul$nb_cell, function(x) tryCatch(t.test(x)[["conf.int"]][2], error=function(e) 0))
    
  }
  
  enh_prop_nb_cell <- list(obs=obs_enh_cell, simul=simul_enh_cell)
  
  ################################################# Save RData ################################################# 
  
  save(enh_prop, enh_prop_nb_cell, enh_prop_dist, file = paste(pathFigures, "RData/data.enhancer.coverage.", ref_sp, ".Rdata", sep=""))

  ############################################################################################################## 

}

############################################################################################################## 
