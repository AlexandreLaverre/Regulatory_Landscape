#########################################################################################################################

library(Hmisc)

options(stringsAsFactors = FALSE)

source("parameters.R") ##  are defined based on the user name

##################################  Correlation gene expression & enhancers activity ####################################

for(ref_sp in c("human", "mouse")){
  enhancers = enhancer.datasets[[ref_sp]]
  
  obs_correl_activity_dist <- list()
  simul_correl_activity_dist <- list()
  
  for (enh in enhancers){
    if(file.exists(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_real_data.txt", sep="/")) & file.exists(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_simulated_data.txt", sep="/"))){
      
      obs <- read.table(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_real_data.txt", sep="/"), h=T, sep="\t")
      simul <- read.table(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_simulated_data.txt", sep="/"), h=T, sep="\t")
      
      # according to distance
    
      obs$dist_class <-cut(obs$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)
      simul$dist_class <- cut(simul$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)
      
      obs_correl_activity_dist[[paste0(enh)]] <- tapply(obs$SpearmanCorrelation, obs$dist_class, function(x) mean(x, na.rm=T))
      obs_correl_activity_dist[[paste0(enh, "_conflow")]] <-tapply(obs$SpearmanCorrelation, obs$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][1])
      obs_correl_activity_dist[[paste0(enh, "_confup")]] <-tapply(obs$SpearmanCorrelation, obs$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][2])
      
      simul_correl_activity_dist[[paste0(enh)]] <- tapply(simul$SpearmanCorrelation, simul$dist_class, function(x) mean(x, na.rm=T))
      simul_correl_activity_dist[[paste0(enh, "_conflow")]] <-tapply(simul$SpearmanCorrelation, simul$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][1])
      simul_correl_activity_dist[[paste0(enh, "_confup")]] <-tapply(simul$SpearmanCorrelation, simul$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][2])
    }
  }
  
  correl_activity <- list(obs=obs_correl_activity_dist, simul=simul_correl_activity_dist)
  
################################################# Save RData ################################################# 
  
  save(correl_activity, file = paste(pathFigures, "RData/data.promoter.enhancer.correlation.", ref_sp, ".Rdata", sep=""))
 
 #########################################################################################################################
}

#########################################################################################################################
