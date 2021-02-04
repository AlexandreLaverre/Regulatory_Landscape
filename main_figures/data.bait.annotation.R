#################################################################################

library(data.table)

source("parameters.R")

#################################################################################

genomes=c("mm10", "hg38")
names(genomes)=c("mouse", "human")

#################################################################################

bait.info=list()
baited.genes=list()

for(sp in c("human", "mouse")){
  baits <- fread(paste(pathFinalData, "SupplementaryDataset1/", sp, "/bait_coords_",genomes[sp],".txt", sep=""), header=T, stringsAsFactors=F, sep="\t")
  class(baits) <- "data.frame"
  
  rownames(baits)=baits$ID
  bait.info[[sp]]=baits
  
  this.genes=unique(unlist(lapply(baits$geneID, function(x) unlist(strsplit(x, split=",")))))
  baited.genes[[sp]]=this.genes
}

#################################################################################

save(list=c("bait.info", "baited.genes"), file=paste(pathFigures, "RData/data.bait.annotation.RData",sep=""))

#################################################################################


