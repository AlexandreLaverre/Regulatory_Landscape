###########################################################################

library(data.table)

source("parameters.R")

pathContacts=paste(pathFinalData, "SupplementaryDataset4/", sep="")

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

###########################################################################

gene.enhancer.contacts=list()

for(sp in c("human", "mouse")){

  info=sampleinfo[[sp]]
  rownames(info)=info$Sample.ID
  
  samples=info$Sample.ID 
  celltypes=info$Broad.cell.type.or.tissue
  names(celltypes)=samples
  
  gene.enhancer.contacts[[sp]]=list()
  
  for(enh in enhancer.datasets[[sp]]){

    real <- fread(paste(pathContacts, sp, "/", enh, "/gene_enhancer_contacts_original_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    class(real) <- "data.frame"
    
    real$nb_cell <- apply(real[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
   
    sim<-fread(paste(pathContacts, sp, "/", enh, "/gene_enhancer_contacts_simulated_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    class(sim) <- "data.frame"
    
    sim$nb_cell <- apply(sim[,samples], 1, function(x) length(unique(celltypes[which(!is.na(x))])))
    
    ## select interactions in the accepted distance range
    real=real[which(real$dist>=minDistance & real$dist<=maxDistance),]
    sim=sim[which(sim$dist>=minDistance & sim$dist<=maxDistance),]

    ## save data
    
    gene.enhancer.contacts[[sp]][[enh]]=list("real"=real, "simulated"=sim)
  }
}

###########################################################################

save(list=c("gene.enhancer.contacts"), file=paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))

###########################################################################


