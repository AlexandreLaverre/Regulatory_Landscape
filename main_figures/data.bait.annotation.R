#################################################################################

library(data.table)

source("parameters.R")

#################################################################################

genomes=c("mm10", "hg38")
names(genomes)=c("mouse", "human")

#################################################################################

bait.info=list()

for(sp in c("human", "mouse")){
  baits <- fread(paste(pathFinalData, "SupplementaryDataset1/", sp, "/bait_coords_",genomes[sp],".txt", sep=""), header=T, stringsAsFactors=F, sep="\t")

  rownames(baits)=baits$ID

  bait.info[[sp]]=baits
}

#################################################################################

save(list=c("bait.info"), file=paste(pathFigures, "RData/data.bait.annotation.RData",sep=""))

#################################################################################


