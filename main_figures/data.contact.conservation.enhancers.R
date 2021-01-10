#######################################################################################

library(data.table)

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7", sep="")

#######################################################################################

contact.conservation=list()

for(ref in c("human", "mouse")){

  tg=setdiff(c("human", "mouse"), ref)

  contact.conservation[[paste(ref, "2", tg, sep="")]]=list()

  for(enh in enhancer.datasets[[ref]]){
    message(enh)
    obsobs=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_original2", tg,"_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    simsim=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_simulated2", tg,"_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    class(obsobs)="data.frame"
    class(simsim)="data.frame"
    
    ## take only orthologous genes presents in both species datasets
    obsobs=obsobs[which(obsobs$target_data == "TRUE"),]
    simsim=simsim[which(simsim$target_data == "TRUE"),]
    
    ## filtered enhancers - single BLAT match, within accepted distance range
    obsobs=obsobs[which(obsobs$BLAT_match == 1 & obsobs$origin_dist >= minDistance & obsobs$origin_dist <= maxDistance),]
    simsim=simsim[which(simsim$BLAT_match == 1 & simsim$origin_dist >= minDistance & simsim$origin_dist <= maxDistance),]

    ## take only well conserved enhancers
    align.threshold <- 0.4 

    obsobs=obsobs[which(obsobs$align_score>=align.threshold),]
    simsim=simsim[which(simsim$align_score>=align.threshold),]
     
    ## save data
    
    contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]]=list("obsobs"=obsobs, "simsim"=simsim) 
    
  }
}

#######################################################################################

save(list=c("contact.conservation"), file=paste(pathFigures, "RData/data.contact.conservation.enhancers.RData",sep=""))

#######################################################################################

