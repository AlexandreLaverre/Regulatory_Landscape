#######################################################################################

library(data.table)

options(stringsAsFactors = FALSE)

source("parameters.R") ## paths are defined based on the user name

path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

#######################################################################################

types=c("pcidentical", "pcungapped")

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep="")) 
load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep="")) ## enhancers are already filtered for duplication levels, repeat proportion etc
load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep="")) ## fragments are already filtered for duplication levels, repeat proprtion etc

#######################################################################################

for(ref_sp in c("human", "mouse")){

  ## bait-fragment contacts, only in cis, unbaited, within a pre-defined distance range
  
  frag.contacts.obs <- observed.contacts[[ref_sp]]
  frag.contacts.sim <- simulated.contacts[[ref_sp]]

  frag.mediandist.obs <- tapply(frag.contacts.obs$distance, as.factor(frag.contacts.obs$id_frag), median)
  frag.mediandist.sim <- tapply(frag.contacts.sim$distance, as.factor(frag.contacts.sim$id_frag), median)
  
  for(type in types){
    
    print(ref_sp)
    
    if(ref_sp == "human"){
      target_sp = "mouse"
      closest_sp="macaque"
    } else{
      target_sp = "human"
      closest_sp="rat"
    }
   
    species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
    
    ## statistics for the contacted restriction fragments
    
    obs <- fragment.statistics[[ref_sp]][["original"]]
    simul <- fragment.statistics[[ref_sp]][["simulated"]]

    ## select fragments that are in previously filtered contacts (within correct distance range)

    obs <- obs[which(obs$ID%in%frag.contacts.obs$id_frag),]
    simul <- simul[which(simul$ID%in%frag.contacts.sim$id_frag),]
    
    ## load alignment score vs all species, for restriction fragments

    frag_align=data.frame("ID"=unique(c(frag.contacts.obs$id_frag, frag.contacts.sim$id_frag)), stringsAsFactors=F)
    rownames(frag_align)=frag_align$ID

    for(tg in species){
      print(tg)
      load(paste(pathFigures, "RData/data.sequence.conservation.fragments.",ref_sp,"2", tg,".RData", sep=""))

      this.cons=get(type)

      frag_align[,tg]=this.cons[rownames(frag_align)]
      
      rm(list=c("pcidentical", "pcungapped"))
    }
    
    frag_align_obs <- frag_align[which(frag_align$ID %in% obs$ID), c("ID", species)]
    frag_align_simul <- frag_align[which(frag_align$ID %in% simul$ID), c("ID", species) ]
    
    ## we compute median distance on filtered contacts

    frag_align_obs$median_dist <- frag.mediandist.obs[frag_align_obs$ID]
    frag_align_simul$median_dist <- frag.mediandist.sim[frag_align_simul$ID]
    
    frag_align_obs$dist_class <- cut(frag_align_obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    frag_align_simul$dist_class <-cut(frag_align_simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    
    ## alignment for enhancers
    
    enhancers=enhancer.datasets[[ref_sp]]
        
    list_align_enh <- list()
    
    for (enh in enhancers){
      print(enh)

      ## filtered contacts for this enhancer set

      enh.contacts.obs <- gene.enhancer.contacts[[ref_sp]][[enh]][["real"]]
      enh.contacts.sim <- gene.enhancer.contacts[[ref_sp]][[enh]][["simulated"]]

      ## median distance for filtered contacts
      
      enh.mediandist.obs <- tapply(enh.contacts.obs$dist, as.factor(enh.contacts.obs$enhancer), median)
      enh.mediandist.sim <- tapply(enh.contacts.sim$dist, as.factor(enh.contacts.sim$enhancer), median)

      enh_align=data.frame("ID"=unique(c(enh.contacts.obs$enh, enh.contacts.sim$enh)), stringsAsFactors=F)
      rownames(enh_align)=enh_align$ID
      
      for(tg in species){
        print(tg)
        load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref_sp,"2", tg,".RData", sep=""))
        
        this.cons=get(type) ## ungapped or identical
        
        enh_align[,tg]=this.cons[rownames(enh_align)]
        
        rm(list=c("pcidentical", "pcungapped"))
      }
              
      ## statistics for enhancers 
      enh_obs_stats <- enhancer.statistics[[ref_sp]][[enh]][["original"]]
      enh_simul_stats <- enhancer.statistics[[ref_sp]][[enh]][["simulated"]]
      
      ## select previously filtered enhancers, in contact within acceptable distances

      enh_obs_stats=enh_obs_stats[which(enh_obs_stats$enh%in%enh.contacts.obs$enhancer),]
      enh_simul_stats=enh_simul_stats[which(enh_simul_stats$enh%in%enh.contacts.sim$enhancer),]
      

      ## select species
      enh_align <- enh_align[, c("ID", species)]
      
      enh_align_obs <- enh_align[which(enh_align$ID %in% enh_obs_stats$enh),]
      enh_align_simul<- enh_align[which(enh_align$ID %in% enh_simul_stats$enh),]
  
      
      ## median distance computed from filtered contacts
      
      enh_align_obs$median_dist <- enh.mediandist.obs[enh_align_obs$ID]
      enh_align_simul$median_dist <- enh.mediandist.sim[enh_align_simul$ID]
      
      enh_align_obs$dist_class <- cut(enh_align_obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      enh_align_simul$dist_class <- cut(enh_align_simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)

      ## store results
      list_align_enh[[enh]] <- list("enh_align_obs"=enh_align_obs, "enh_align_simul"=enh_align_simul)
    }

    ## save results

    save(list=c("frag_align_obs", "frag_align_simul", "list_align_enh"), file=paste(pathFigures, "RData/data.sequence.conservation.stats.",type,".",ref_sp,".RData", sep=""))

  }
}

################################################################################################################################################
