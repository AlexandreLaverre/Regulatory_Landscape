#########################################################################################################################

library(Hmisc)
library(data.table)

options(stringsAsFactors = FALSE)

source("parameters.R") ##  are defined based on the user name

load(paste(pathFigures,  "RData/data.fragment.contacts.RData",sep=""))

##################################  Correlation gene expression & enhancers activity ####################################

for(ref_sp in c("human", "mouse")){

  ## fragment contacts
  
  frag.contact.obs=observed.contacts[[ref_sp]]
  frag.contact.sim=simulated.contacts[[ref_sp]]

  frag.contact.obs$idcontact=paste(frag.contact.obs$id_bait, frag.contact.obs$id_frag, sep="-")
  frag.contact.sim$idcontact=paste(frag.contact.sim$id_bait, frag.contact.sim$id_frag, sep="-")
  
  enhancers = enhancer.datasets[[ref_sp]]
  
  obs_correl_activity_dist <- list()
  simul_correl_activity_dist <- list()
  
  for (enh in enhancers){
    if(file.exists(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_real_data.txt", sep="/")) & file.exists(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_simulated_data.txt", sep="/"))){
      
      obs <- fread(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_real_data.txt", sep="/"), h=T, sep="\t")
      class(obs)<-"data.frame"

      ## keep only previously filtered bait-fragment contacts
      obs$IDBait=unlist(lapply(obs$IDBait, function(x) paste(unlist(strsplit(x, split=",")), collapse=":")))
      obs$IDContactedFragment=unlist(lapply(obs$IDContactedFragment, function(x) paste(unlist(strsplit(x, split=",")),collapse=":")))
      obs$IDContact=paste(obs$IDBait, obs$IDContactedFragment, sep="-")

      obs=obs[which(obs$IDContact%in%frag.contact.obs$idcontact),]
      
      simul <- fread(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_simulated_data.txt", sep="/"), h=T, sep="\t")
      class(simul)<-"data.frame"

      ## keep only previously filtered bait-fragment contacts
      simul$IDBait=unlist(lapply(simul$IDBait, function(x) paste(unlist(strsplit(x, split=",")), collapse=":")))
      simul$IDContactedFragment=unlist(lapply(simul$IDContactedFragment, function(x) paste(unlist(strsplit(x, split=",")), collapse=":")))
      simul$IDContact=paste(simul$IDBait, simul$IDContactedFragment, sep="-")

      simul=simul[which(simul$IDContact%in%frag.contact.sim$idcontact),]

      print(paste(ref_sp, enh, "observed mean corr, all distances", mean(obs$SpearmanCorrelation)))
      print(paste(ref_sp, enh, "simulated mean corr, all distances", mean(simul$SpearmanCorrelation)))
      print(wilcox.test(obs$SpearmanCorrelation, simul$SpearmanCorrelation))
            
      # according to distance
    
      obs$dist_class <-cut(obs$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)
      simul$dist_class <- cut(simul$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)

      obs_correl_activity_dist[[paste0(enh)]] <- tapply(obs$SpearmanCorrelation, obs$dist_class, function(x) mean(x, na.rm=T))
      obs_correl_activity_dist[[paste0(enh, "_conflow")]] <- tapply(obs$SpearmanCorrelation, obs$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][1])
      obs_correl_activity_dist[[paste0(enh, "_confup")]] <- tapply(obs$SpearmanCorrelation, obs$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][2])

      simul_correl_activity_dist[[paste0(enh)]] <- tapply(simul$SpearmanCorrelation, simul$dist_class, function(x) mean(x, na.rm=T))
      simul_correl_activity_dist[[paste0(enh, "_conflow")]] <- tapply(simul$SpearmanCorrelation, simul$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][1])
      simul_correl_activity_dist[[paste0(enh, "_confup")]] <- tapply(simul$SpearmanCorrelation, simul$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][2])
    }
  }
  
  correl_activity <- list(obs=obs_correl_activity_dist, simul=simul_correl_activity_dist)
  
################################################# Save RData ################################################# 
  
  save(correl_activity, file = paste(pathFigures, "RData/data.promoter.enhancer.correlation.", ref_sp, ".RData", sep=""))
 
 #########################################################################################################################
}

#########################################################################################################################
