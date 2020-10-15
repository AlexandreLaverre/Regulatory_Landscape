###########################################################################

source("parameters.R")

pathContacts=paste(pathFinalData, "SupplementaryDataset4/", sep="")

###########################################################################

enhancer.statistics=list()

for(sp in c("human", "mouse")){
  
  enhancer.statistics[[sp]]=list()
  
  for(enh in enhancer.datasets[[sp]]){
    
    real=read.table(paste(pathContacts, sp, "/", enh, "/statistics_contacted_enhancers_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    sim=read.table(paste(pathContacts, sp, "/", enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    
    enhancer.statistics[[sp]][[enh]]=list("real"=real, "simulated"=sim)
  }
}

###########################################################################

save(list=c("enhancer.statistics"), file=paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))

###########################################################################
