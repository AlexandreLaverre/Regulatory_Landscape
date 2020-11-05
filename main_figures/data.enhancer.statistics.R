###########################################################################

library(data.table)

source("parameters.R")

pathStats=paste(pathFinalData, "SupplementaryDataset4/", sep="")

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
    
    ## select enhancers that are not duplicated and with repeat_part < 1%
    
    obs <- obs[which(obs$BLAT_match < 2 & (obs$repeat_bp/obs$length) < 0.01),] 
    sim <- sim[which(sim$BLAT_match < 2 & (sim$repeat_bp/sim$length) < 0.01),]
    
    ## save results
    
    enhancer.statistics[[sp]][[enh]]=list("original"=obs, "simulated"=sim)
  }
}

###########################################################################

save(list=c("enhancer.statistics"), file=paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))

###########################################################################
