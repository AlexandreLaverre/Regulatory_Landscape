#######################################################################################

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7", sep="")

#######################################################################################

contact.conservation=list()

for(ref in c("human", "mouse")){

  tg=setdiff(c("human", "mouse"), ref)

  contact.conservation[[paste(ref, "2", tg, sep="")]]=list()
  
  for(enh in enhancer.datasets[[ref]]){

    obsobs=read.table(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_original2", tg,"_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    obssim=read.table(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_original2", tg,"_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    simobs=read.table(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_simulated2", tg,"_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    simsim=read.table(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_simulated2", tg,"_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]]=list("obsobs"=obsobs, "simsim"=simsim, "obssim"=obssim, "simobs"=simobs)
    
  }
}

#######################################################################################

save(list=c("contact.conservation"), file=paste(pathFigures, "RData/data.contact.conservation.RData",sep=""))

#######################################################################################

