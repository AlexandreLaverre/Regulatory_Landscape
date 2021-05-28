###########################################################################

library(data.table)

source("../main_figures/parameters.R")

pathFragments=paste(pathFinalData, "SupplementaryDataset1/", sep="")

###########################################################################

species=c("human", "mouse")
genomes=c("hg38", "mm10")

names(genomes)=species

###########################################################################

restriction.map=list()

###########################################################################

for(sp in species){
  print(sp)
  
  genome=genomes[sp]
  
  frag=fread(paste(pathFragments, sp, "/frag_coords_",genome,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  class(frag)<-"data.frame"

  restriction.map[[sp]]=frag
}

###########################################################################

save(list=c("restriction.map"), file=paste(pathFigures, "RData/data.restriction.map.RData", sep=""))

###########################################################################
