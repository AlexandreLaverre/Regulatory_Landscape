###########################################################################

source("parameters.R")

pathNeighborEnhancers="../../results/neighbor_enhancers/"

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
  
  for(enh in enhancer.datasets[[sp]]){
    
    print(enh)
    
    real=read.table(paste(pathNeighborEnhancers, sp, "/", enh, "/neighbor_enhancers_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F)
    
    ## select protein-coding genes

    real=real[which(real$GeneID%in%pc.genes),]

    ## select enhancers in our accepted distance range
       
    real=real[which(real$Distance>=minDistance & real$Distance<=maxDistance),]

    ## change column names

    colnames(real)[which(colnames(real)=="EnhancerID")]="enhancer"
    colnames(real)[which(colnames(real)=="GeneID")]="gene"

    ## add chr to enhancer id

    g=grep("^chr", real$enhancer)

    if(length(g)==0){
      real$enhancer=paste("chr", real$enhancer, sep="")
    }
    ## save data
    
    gene.enhancer.contacts[[sp]][[enh]]=list("real"=real)
  }
}

###########################################################################

save(list=c("gene.enhancer.contacts"), file=paste(pathFigures, "RData/data.neighbor.enhancers.RData", sep=""))

###########################################################################


