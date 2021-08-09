#########################################################################################################################

source("parameters.R") ##  are defined based on the user name

library(Hmisc)
library(data.table)
library(bootBCa, lib=pathRlibs)

set.seed(19)

options(stringsAsFactors = FALSE)

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
  neighbors_correl_activity_dist <- list()
  
  for (enh in enhancers){
    if(file.exists(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_real_data.txt", sep="/")) & file.exists(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_simulated_data.txt", sep="/"))){

      print(paste(ref_sp, enh))
      
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

      ## neighbor enhancers

       neighbors <- fread(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_neighbor_promoters_enhancers.txt", sep="/"), h=T, sep="\t")
      class(neighbors)<-"data.frame"

      ## filter for distances

      neighbors=neighbors[which(neighbors$Distance>=minDistance & neighbors$Distance<=maxDistance),]

      ## select the same promoters as for the observed contacts

      neighbors=neighbors[which(neighbors$IDPromoter%in%obs$IDPromoter),]

      print(paste(ref_sp, enh, "observed mean corr, all distances", mean(obs$SpearmanCorrelation)))
      print(paste(ref_sp, enh, "simulated mean corr, all distances", mean(simul$SpearmanCorrelation)))
      print(paste(ref_sp, enh, "neighbors mean corr, all distances", mean(neighbors$SpearmanCorrelation)))
      print(wilcox.test(obs$SpearmanCorrelation, simul$SpearmanCorrelation))
            
      # according to distance
    
      obs$dist_class <-cut(obs$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)
      simul$dist_class <- cut(simul$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)
      neighbors$dist_class <- cut(neighbors$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)

      BC.obs=tapply(obs$SpearmanCorrelation, obs$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
      obs_correl_activity_dist[[enh]]=unlist(lapply(BC.obs, function(x) x[3]))
      obs_correl_activity_dist[[paste0(enh, "_conflow")]]=unlist(lapply(BC.obs, function(x) x[4]))
      obs_correl_activity_dist[[paste0(enh, "_confup")]]=unlist(lapply(BC.obs, function(x) x[5]))
       
      BC.simul=tapply(simul$SpearmanCorrelation, simul$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
      simul_correl_activity_dist[[enh]]=unlist(lapply(BC.simul, function(x) x[3]))
      simul_correl_activity_dist[[paste0(enh, "_conflow")]]=unlist(lapply(BC.simul, function(x) x[4]))
      simul_correl_activity_dist[[paste0(enh, "_confup")]]=unlist(lapply(BC.simul, function(x) x[5]))

      BC.neighbors=tapply(neighbors$SpearmanCorrelation, neighbors$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
      neighbors_correl_activity_dist[[enh]]=unlist(lapply(BC.neighbors, function(x) x[3]))
      neighbors_correl_activity_dist[[paste0(enh, "_conflow")]]=unlist(lapply(BC.neighbors, function(x) x[4]))
      neighbors_correl_activity_dist[[paste0(enh, "_confup")]]=unlist(lapply(BC.neighbors, function(x) x[5]))
    }
  }
  
  correl_activity <- list(obs=obs_correl_activity_dist, simul=simul_correl_activity_dist, neighbors=neighbors_correl_activity_dist)
  
################################################# Save RData ################################################# 
  
  save(correl_activity, file = paste(pathFigures, "RData/data.promoter.enhancer.correlation.", ref_sp, ".RData", sep=""))
 
 #########################################################################################################################
}

#########################################################################################################################
