##################################################################

source("parameters.R")

pathAnnot=paste(pathFinalData, "SupplementaryDataset3/genes/",sep="")

##################################################################

gene.annot=list()

for(sp in c("human", "mouse")){
  this.annot=read.table(paste(pathAnnot, sp, "_genes_Ensembl94.txt",sep=""), h=T, sep="\t", quote=F)
  rownames(this.annot)=this.annot$GeneID

  gene.annot[[sp]]=this.annot
}

##################################################################

save(list=c("gene.annot"), file=paste(pathFigures, "RData/data.gene.annotations.RData",sep=""))

##################################################################
