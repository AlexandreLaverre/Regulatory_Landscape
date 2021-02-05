#########################################################################################################################

library(Hmisc)
library(data.table)

options(stringsAsFactors = FALSE)

source("parameters.R") ##  are defined based on the user name

#########################################################################################################################

load(paste(pathFigures,  "RData/data.fragment.contacts.RData",sep=""))

minTPM=1

#########################################################################################################################

##################################  number of samples with shared expression ############################################

for(ref_sp in c("human", "mouse")){

  ## fragment contacts
  
  frag.contact.obs=observed.contacts[[ref_sp]]
  frag.contact.sim=simulated.contacts[[ref_sp]]

  frag.contact.obs$idcontact=paste(frag.contact.obs$id_bait, frag.contact.obs$id_frag, sep="-")
  frag.contact.sim$idcontact=paste(frag.contact.sim$id_bait, frag.contact.sim$id_frag, sep="-")
  
  enhancers = enhancer.datasets[[ref_sp]]
  
  obs_shared_activity_dist <- list()
  simul_shared_activity_dist <- list()
  
  for (enh in enhancers){
    pathObs=paste(pathFinalData, "SupplementaryDataset8/", ref_sp, "/", enh, "/shared_activity_minTPM",minTPM,"_real_data.txt", sep="")
    pathSim=paste(pathFinalData, "SupplementaryDataset8/", ref_sp, "/", enh, "/shared_activity_minTPM",minTPM,"_simulated_data.txt", sep="")
    
    if(file.exists(pathObs) & file.exists(pathSim)){
      
      obs <- fread(pathObs, h=T, sep="\t")
      class(obs)<-"data.frame"

      ## keep only previously filtered bait-fragment contacts
      obs$IDBait=unlist(lapply(obs$IDBait, function(x) paste(unlist(strsplit(x, split=",")), collapse=":")))
      obs$IDContactedFragment=unlist(lapply(obs$IDContactedFragment, function(x) paste(unlist(strsplit(x, split=",")),collapse=":")))
      obs$IDContact=paste(obs$IDBait, obs$IDContactedFragment, sep="-")

      obs=obs[which(obs$IDContact%in%frag.contact.obs$idcontact),]
      
      simul <- fread(pathSim, h=T, sep="\t")
      class(simul)<-"data.frame"

      ## keep only previously filtered bait-fragment contacts
      simul$IDBait=unlist(lapply(simul$IDBait, function(x) paste(unlist(strsplit(x, split=",")), collapse=":")))
      simul$IDContactedFragment=unlist(lapply(simul$IDContactedFragment, function(x) paste(unlist(strsplit(x, split=",")), collapse=":")))
      simul$IDContact=paste(simul$IDBait, simul$IDContactedFragment, sep="-")

      simul=simul[which(simul$IDContact%in%frag.contact.sim$idcontact),]

      ## fraction of samples with shared activity

      obs$PropShared=obs$NbExpBoth/obs$NbExpPromoter
      simul$PropShared=simul$NbExpBoth/simul$NbExpPromoter
            
      # according to distance
    
      obs$dist_class <-cut(obs$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)
      simul$dist_class <- cut(simul$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)

      obs_shared_activity_dist[[enh]] <- tapply(obs$PropShared, obs$dist_class, function(x) mean(x, na.rm=T))
      obs_shared_activity_dist[[paste0(enh, "_conflow")]] <- tapply(obs$PropShared, obs$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][1])
      obs_shared_activity_dist[[paste0(enh, "_confup")]] <- tapply(obs$PropShared, obs$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][2])

      simul_shared_activity_dist[[enh]] <- tapply(simul$PropShared, simul$dist_class, function(x) mean(x, na.rm=T))
      simul_shared_activity_dist[[paste0(enh, "_conflow")]] <- tapply(simul$PropShared, simul$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][1])
      simul_shared_activity_dist[[paste0(enh, "_confup")]] <- tapply(simul$PropShared, simul$dist_class, function(x) t.test(x, na.rm=T)[["conf.int"]][2])

    
  }
  
  shared_activity <- list(obs=obs_shared_activity_dist, simul=simul_shared_activity_dist)
  
################################################# Save RData ################################################# 
  
  save(shared_activity, file = paste(pathFigures, "RData/data.promoter.enhancer.shared.activity.", ref_sp, ".RData", sep=""))
 
 #########################################################################################################################
}

#########################################################################################################################
