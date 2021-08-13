###########################################################################

source("../main_figures/parameters.R")

pathNeighborEnhancers=paste(pathFinalData, "SupplementaryDataset4/",sep="")

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))

###########################################################################

gene.enhancer.contacts=list()

for(sp in c("human", "mouse")){
  print(sp)

  ## gene annotations 

  annot=gene.annot[[sp]]
  pc.genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
  
  ## gene enhancer contacts
  
  gene.enhancer.contacts[[sp]]=list()

  regions=read.table(paste(pathNeighborEnhancers, sp, "/predicted_regulatory_regions_neighbor_TSS.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

  gene.enhancer.contacts[[sp]][["regions"]]=regions
  
  for(enh in enhancer.datasets[[sp]]){
    
    print(enh)
    
    real=read.table(paste(pathNeighborEnhancers, sp, "/", enh, "/predicted_enhancers_neighbor_TSS.txt", sep=""), h=T, stringsAsFactors=F)
    
    ## select protein-coding genes

    real=real[which(real$GeneID%in%pc.genes),] 
       
    ## save data
    
    gene.enhancer.contacts[[sp]][[enh]]=list("real"=real)
  }
}

###########################################################################

save(list=c("gene.enhancer.contacts"), file=paste(pathFigures, "RData/data.neighbor.enhancers.RData", sep=""))

###########################################################################


