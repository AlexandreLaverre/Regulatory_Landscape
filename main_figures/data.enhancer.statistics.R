###########################################################################

library(data.table)

source("parameters.R")

pathStats=paste(pathFinalData, "SupplementaryDataset4/", sep="")

load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep="")) ## fragments are already filtered for duplication levels, repeat proprtion etc

###########################################################################

enhancer.statistics=list()

for(sp in c("human", "mouse")){
  
  enhancer.statistics[[sp]]=list()
  
  for(enh in enhancer.datasets[[sp]]){
    
    obs=fread(paste(pathStats, sp, "/", enh, "/statistics_contacted_enhancers_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    sim=fread(paste(pathStats, sp, "/", enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

    class(obs)<-"data.frame"
    class(sim)<-"data.frame"

    ## create enhancer identifier
    
    obs$enh <-  do.call(paste,c(obs[c("chr","start","end")],sep=":"))
    sim$enh <-  do.call(paste,c(sim[c("chr","start","end")],sep=":"))

    ## rownames
    rownames(obs) <- obs[,"enh"]
    rownames(sim) <- sim[,"enh"]
    
    print(paste(nrow(obs), "observed contacted enhancers"))
    print(paste(nrow(sim), "simulated contacted enhancers"))
    
    ## select enhancers depending on the number of BLAT hits
    
    obs <- obs[which(obs$BLAT_match > minBLAT & obs$BLAT_match < maxBLAT),]
    sim <- sim[which(sim$BLAT_match > minBLAT & sim$BLAT_match < maxBLAT),]
    
    ## select enhancers within previously filtered fragments
    enh.to.frag = fread(paste(pathStats, sp, "/", enh, "/enhancers_to_restriction_fragment_IDs.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    
    frag.obs <- fragment.statistics[[sp]][["original"]]
    frag.sim <- fragment.statistics[[sp]][["simulated"]]
    
    enh.obs.to.keep <- enh.to.frag[which(enh.to.frag$ID.fragment %in% rownames(frag.obs)),]
    enh.sim.to.keep <- enh.to.frag[which(enh.to.frag$ID.fragment %in% rownames(frag.sim)),]
    
    obs <- obs[which(obs$enh %in% enh.obs.to.keep$ID.enhancer),]
    sim <- sim[which(sim$enh %in% enh.sim.to.keep$ID.enhancer),]
    
  
    print(paste(nrow(obs), "observed contacted enhancers after filtering"))
    print(paste(nrow(sim), "simulated contacted enhancers after filtering"))
    
    ## save results
    
    enhancer.statistics[[sp]][[enh]]=list("original"=obs, "simulated"=sim)
  }
}

###########################################################################

save(list=c("enhancer.statistics"), file=paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))

###########################################################################
