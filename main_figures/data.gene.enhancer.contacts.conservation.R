###########################################################################

source("parameters.R")

pathConservation=paste(pathFinalData, "SupplementaryDataset7/", sep="")

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

###########################################################################

contacts.conservation=list()

for(sp in c("human", "mouse")){
  
  info=sampleinfo[[sp]]
  rownames(info)=info$Sample.ID
  
  samples=info$Sample.ID 
  celltypes=info$Broad.cell.type.or.tissue
  names(celltypes)=samples
  
  if (sp == "human"){target_sp="mouse"}else{target_sp="human"}
  contacts.conservation[[sp]]=list()
  
  for(enh in enhancer.datasets[[sp]]){

    real=read.table(paste(pathConservation, sp, "/contact_conservation/", enh, "/", sp, "_original2", target_sp, "_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    real$nb_cell <- apply(real[,samples],1, function(x) length(unique(celltypes[which(x>0)])))
    
    sim=read.table(paste(pathConservation, sp, "/contact_conservation/", enh, "/", sp, "_simulated2", target_sp, "_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    sim$nb_cell <- apply(sim[,samples], 1, function(x) length(unique(celltypes[which(x>0)])))
    
    contacts.conservation[[sp]][[enh]]=list("real"=real, "simulated"=sim)
  }
}

###########################################################################

save(list=c("contacts.conservation"), file=paste(pathFigures, "RData/data.gene.enhancer.contacts.conservation.RData", sep=""))

###########################################################################
