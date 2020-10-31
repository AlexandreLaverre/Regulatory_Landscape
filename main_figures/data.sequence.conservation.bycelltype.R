################################################################################################################################################

library(data.table)

options(stringsAsFactors = FALSE)

source("parameters.R") ## paths are defined based on the user name

outnames=c("pcidentity", "pcungapped")
names(outnames)=c("IdenticalSequence", "UngappedAlignment")

## load interaction data

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep="")) 
load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))

################################################################################################################################################

for(ref_sp in c("human", "mouse")){
  
  ## bait-fragment contacts, only in cis, unbaited, max distance 2Mb, min distance 25kb
  
  frag.contacts.obs <- observed.contacts[[ref_sp]]
  frag.contacts.sim <- simulated.contacts[[ref_sp]]

  ## sample info

  info=sampleinfo[[ref_sp]]

  cells=unique(info[,"Broad.cell.type.or.tissue"])
  cell.samples <- tapply(info$Sample.ID, as.factor(info$Broad.cell.type.or.tissue), function(x) x)

  for(c in cells){
    frag.contacts.obs[,c]=apply(as.matrix(frag.contacts.obs[,cell.samples[[c]]]), 1, function(x) any(!is.na(x)))
    frag.contacts.sim[,c]=apply(as.matrix(frag.contacts.sim[,cell.samples[[c]]]), 1, function(x) any(!is.na(x)))
  }

  ## cell types contacted by each restriction fragment

  frag.cells.obs <- apply(frag.contacts.obs[,cells], 2, function(x) tapply(x, as.factor(frag.contacts.obs$id_frag), any))
  frag.cells.sim <- apply(frag.contacts.sim[,cells], 2, function(x) tapply(x, as.factor(frag.contacts.sim$id_frag), any))
  
  ## alignment statistics
  
  for(type in c("IdenticalSequence", "UngappedAlignment")){
    
    print(ref_sp)
    
    if(ref_sp == "human"){
      target_sp = "mouse"
      closest_sp="macaque"
    } else{
      target_sp = "human"
      closest_sp="rat"
    }

    enhancers=enhancer.datasets[[ref_sp]]
    
    species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
    
    path_evol <- paste(pathFinalData, "SupplementaryDataset7/", ref_sp, "/", sep="")
    path_annot <- paste(pathFinalData, "SupplementaryDataset4/", ref_sp, "/", sep="")
    
### statistics for the contacted restriction fragments

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

    ## select previously filtered contacts

    obs <- obs[which(obs$ID%in%frag.contacts.obs$id_frag),]
    simul <- simul[which(simul$ID%in%frag.contacts.sim$id_frag),]
    
    ## alignment score vs all species, for restriction fragments 
    frag_align <- read.table(paste(path_evol,"/sequence_conservation/restriction_fragments/AlignmentStatistics_Excluding_Exons_",type,"_AllSpecies.txt", sep=""), header=T)

    ## replace NA values with 0
    for(sp in species){
      frag_align[which(is.na(frag_align[,sp])),sp]=0
    }
    
    frag_align_obs <- frag_align[which(frag_align$ID %in% obs$ID), c("ID", species)]
    frag_align_simul <- frag_align[which(frag_align$ID %in% simul$ID), c("ID", species) ]
    
    frag_align_obs$median_dist <- obs[frag_align_obs$ID,"median_dist"]
    frag_align_simul$median_dist <- simul[frag_align_simul$ID,"median_dist"]

    frag_align_obs$dist_class <- cut(frag_align_obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    frag_align_simul$dist_class <-cut(frag_align_simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)

    ## add cell type info

    frag_align_obs=cbind(frag_align_obs, frag.cells.obs[frag_align_obs$ID,])
    frag_align_simul=cbind(frag_align_simul, frag.cells.sim[frag_align_simul$ID,])
    
    ## alignment for enhancers
    
    list_align_enh <- list()
    
    for (enh in enhancers){
      print(enh)

      ## filtered contacts for this enhancer set

      enh.contacts.obs <- gene.enhancer.contacts[[ref_sp]][[enh]][["real"]]
      enh.contacts.sim <- gene.enhancer.contacts[[ref_sp]][[enh]][["simulated"]]

      ## cell types
      for(c in cells){
        enh.contacts.obs[,c]=apply(as.matrix(enh.contacts.obs[,cell.samples[[c]]]), 1, function(x) any(!is.nan(x)))
        enh.contacts.sim[,c]=apply(as.matrix(enh.contacts.sim[,cell.samples[[c]]]), 1, function(x) any(!is.nan(x)))
      }
      
      ## cell types contacted by each restriction fragment
      
      enh.cells.obs <- apply(enh.contacts.obs[,cells], 2, function(x) tapply(x, as.factor(enh.contacts.obs$enhancer), any))
      enh.cells.sim <- apply(enh.contacts.sim[,cells], 2, function(x) tapply(x, as.factor(enh.contacts.sim$enhancer), any))

      ## alignment scores
      
      enh_align <- read.table(paste(path_evol,"/sequence_conservation/enhancers/", enh, "/AlignmentStatistics_Excluding_Exons_",type,"_AllSpecies.txt", sep=""), header=T)

      ## replace NA values with 0
      for(sp in species){
        enh_align[which(is.na(enh_align[,sp])),sp]=0
      }

      ## statistics for enhancers 
      
      enh_obs_stats <- read.table(paste(path_annot, enh, "/statistics_contacted_enhancers_original.txt", sep=""), header=T)
      enh_simul_stats <- read.table(paste(path_annot, enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), header=T)

      enh_obs_stats$enh <-  do.call(paste,c(enh_obs_stats[c("chr","start","end")],sep=":"))
      enh_simul_stats$enh <-  do.call(paste,c(enh_simul_stats[c("chr","start","end")],sep=":"))

      rownames(enh_obs_stats) <- enh_obs_stats[,"enh"]
      rownames(enh_simul_stats) <- enh_simul_stats[,"enh"]

      ## select previously filtered enhancers

      enh_obs_stats=enh_obs_stats[which(enh_obs_stats$enh%in%enh.contacts.obs$enhancer),]
      enh_simul_stats=enh_simul_stats[which(enh_simul_stats$enh%in%enh.contacts.sim$enhancer),]
      
      ## distance class
      enh_obs_stats$dist_class <- cut(enh_obs_stats$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      enh_simul_stats$dist_class <-cut(enh_simul_stats$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      
      
      ## Select unduplicated and with repeat_part < 1%
      enh_obs_stats <- enh_obs_stats[which(enh_obs_stats$BLAT_match < 2 & (enh_obs_stats$repeat_bp/enh_obs_stats$length) < 0.01),] 
      enh_simul_stats <- enh_simul_stats[which(enh_simul_stats$BLAT_match < 2 & (enh_simul_stats$repeat_bp/enh_simul_stats$length) < 0.01),]
      
      ## select species
      enh_align <- enh_align[, c("ID", species)]
      
      enh_align_obs <- enh_align[which(enh_align$ID %in% enh_obs_stats$enh),]
      enh_align_simul<- enh_align[which(enh_align$ID %in% enh_simul_stats$enh),]

      enh_align_obs$median_dist <- enh_obs_stats[enh_align_obs$ID, "median_dist"]
      enh_align_obs$dist_class <- enh_obs_stats[enh_align_obs$ID, "dist_class"]

      enh_align_simul$median_dist <- enh_simul_stats[enh_align_simul$ID, "median_dist"]
      enh_align_simul$dist_class <- enh_simul_stats[enh_align_simul$ID, "dist_class"]

      ## add cell type info
      
      enh_align_obs=cbind(enh_align_obs, enh.cells.obs[enh_align_obs$ID,])
      enh_align_simul=cbind(enh_align_simul, enh.cells.sim[enh_align_simul$ID,])

      ## store results
      list_align_enh[[enh]] <- list("enh_align_obs"=enh_align_obs, "enh_align_simul"=enh_align_simul)
    }

    ## save results

    save(list=c("frag_align_obs", "frag_align_simul", "list_align_enh"), file=paste(pathFigures, "RData/data.sequence.conservation.celltypes.",outnames[type],".",ref_sp,".Rdata", sep=""))

  }
}

################################################################################################################################################
