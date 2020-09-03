#######################################################################################

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7", sep="")

#######################################################################################

contact.conservation=list()

for(ref in c("human", "mouse")){

  tg=setdiff(c("human", "mouse"), ref)

  contact.conservation[[paste(ref, "2", tg, sep="")]]=list()
  
  for(enh in enhancer.datasets[[ref]]){

    obs=read.table(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "2", tg,"_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    sim=read.table(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "2", tg,"_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]]=list("observed"=obs, "simulated"=sim)
    
  }
}

#######################################################################################

save(list=c("contact.conservation"), file=paste(pathFigures, "RData/data.contact.conservation.RData",sep=""))

#######################################################################################

