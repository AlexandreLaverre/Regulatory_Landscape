#######################################################################################

library(data.table)

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7", sep="")


#######################################################################################

load(paste(pathFigures,"RData/data.enhancer.statistics.RData", sep=""))

#######################################################################################

contact.conservation=list()

for(ref in c("human", "mouse")){

  tg=setdiff(c("human", "mouse"), ref)

  contact.conservation[[paste(ref, "2", tg, sep="")]]=list()

  for(enh in enhancer.datasets[[ref]]){

    obsobs=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_original2", tg,"_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    obssim=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_original2", tg,"_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    simobs=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_simulated2", tg,"_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    simsim=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_simulated2", tg,"_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    class(obsobs)="data.frame"
    class(obssim)="data.frame"
    class(simobs)="data.frame"
    class(simsim)="data.frame"
    
    ## take only orthologous genes presents in both species datasets
    obsobs=obsobs[which(obsobs$target_gene != "NA"),]
    obssim=obssim[which(obssim$target_gene != "NA"),]
    simobs=simobs[which(simobs$target_gene != "NA"),]
    simsim=simsim[which(simsim$target_gene != "NA"),]
    
    ## filtered enhancers
    enh.obs=enhancer.statistics[[ref]][[enh]][["original"]]
    enh.sim=enhancer.statistics[[ref]][[enh]][["simulated"]]

    ## select previously filtered enhancers
    obsobs=obsobs[which(obsobs$origin_enh%in%enh.obs$enh),]
    obssim=obssim[which(obssim$origin_enh%in%enh.obs$enh),]
    simobs=simobs[which(simobs$origin_enh%in%enh.sim$enh),]
    simsim=simsim[which(simsim$origin_enh%in%enh.sim$enh),]

    ## take only well conserved and unduplicated enhancers
    ## threshold alignment score: 5% quantile, observed values
    
    align.threshold <- 0.4 #quantile(obsobs$align_score, p=0.05, na.rm=T)

    obsobs=obsobs[which(obsobs$align_score>=align.threshold & obsobs$BLAT_match < 2),]
    obssim=obssim[which(obssim$align_score>=align.threshold & obssim$BLAT_match < 2),]
    simobs=simobs[which(simobs$align_score>=align.threshold & simobs$BLAT_match < 2),]
    simsim=simsim[which(simsim$align_score>=align.threshold & simsim$BLAT_match < 2),]
   
    
    ## save data
    
    contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]]=list("obsobs"=obsobs, "simsim"=simsim, "obssim"=obssim, "simobs"=simobs)
    
  }
}

#######################################################################################

save(list=c("contact.conservation"), file=paste(pathFigures, "RData/data.contact.conservation.enhancers.RData",sep=""))

#######################################################################################

