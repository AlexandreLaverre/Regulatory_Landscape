###########################################################################

source("parameters.R")

pathContacts=paste(pathFinalData, "SupplementaryDataset4/", sep="")

###########################################################################

gene.enhancer.contacts=list()

for(sp in c("human", "mouse")){

  gene.enhancer.contacts[[sp]]=list()
  
  for(enh in enhancer.datasets[[sp]]){

    real=read.table(paste(pathContacts, sp, "/", enh, "/gene_enhancer_contacts_original_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

    sim=read.table(paste(pathContacts, sp, "/", enh, "/gene_enhancer_contacts_simulated_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

    gene.enhancer.contacts[[sp]][[enh]]=list("real"=real, "simulated"=sim)
  }
}

###########################################################################

save(list=c("gene.enhancer.contacts"), file=paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))

###########################################################################


