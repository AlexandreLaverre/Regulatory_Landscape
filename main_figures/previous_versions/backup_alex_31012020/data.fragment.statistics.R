###########################################################################

library(data.table)

source("parameters.R")

pathStats=paste(pathFinalData, "SupplementaryDataset5/", sep="")

###########################################################################

fragment.statistics=list()

for(sp in c("human", "mouse")){
  
  obs=fread(paste(pathStats, sp, "/statistics_contacted_sequence_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  sim=fread(paste(pathStats, sp, "/statistics_contacted_sequence_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  
  class(obs)<-"data.frame"
  class(sim)<-"data.frame"
  
  ## create fragment identifier
  
  obs$ID <-  do.call(paste,c(obs[c("chr","start","end")],sep=":"))
  sim$ID <-  do.call(paste,c(sim[c("chr","start","end")],sep=":"))
  
  ## rownames
  rownames(obs) <- obs$ID
  rownames(sim) <- sim$ID

   ## unbaited only
  obs <- obs[which(obs$baited == "unbaited"),]
  sim <- sim[which(sim$baited == "unbaited"),]

  ## select fragments that are not duplicated 
  
  obs <- obs[which(obs$BLAT_match < 2),] 
  sim <- sim[which(sim$BLAT_match < 2),]
  
  ## save results
  
  fragment.statistics[[sp]]=list("original"=obs, "simulated"=sim)
}

###########################################################################

save(list=c("fragment.statistics"), file=paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))

###########################################################################
