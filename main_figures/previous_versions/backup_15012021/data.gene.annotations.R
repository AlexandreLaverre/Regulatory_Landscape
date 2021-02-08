##################################################################

library(data.table)

source("parameters.R")

pathAnnot=paste(pathFinalData, "SupplementaryDataset3/genes/",sep="")

##################################################################

gene.annot=list()

for(sp in c("human", "mouse")){
  this.annot=fread(paste(pathAnnot, sp, "_genes_Ensembl94.txt",sep=""), h=T, sep="\t", quote="\"")
  rownames(this.annot)=this.annot$GeneID

  class(this.annot)="data.frame"

  gene.annot[[sp]]=this.annot
}

##################################################################

save(list=c("gene.annot"), file=paste(pathFigures, "RData/data.gene.annotations.RData",sep=""))

##################################################################
