########################################################################

source("parameters.R")

########################################################################

load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

########################################################################

library(ade4)

########################################################################

data.AFC=list()

for(sp in c("human")){
  info=sampleinfo[[sp]]
  samples=info$Sample.ID

  data_AFC <- observed.contacts[[sp]]
  data_AFC <- data_AFC[, samples] 
  data_AFC [!is.na(data_AFC)] <- 1
  data_AFC [is.na(data_AFC)] <- 0
  data_AFC <- data.frame(t(data_AFC)) # row = samples

  this.AFC <- dudi.coa(data_AFC, scannf=F, nf=3)

  dist.AFC=dist.dudi(this.AFC) 
  hclust.AFC=hclust(dist.AFC, "ward")
  
  data.AFC[[sp]]=list("AFC"=this.AFC, "hclust.AFC"=hclust.AFC)
}

########################################################################

save(data.AFC, file=paste(pathFigures, "RData/data.AFC.RData", sep=""))

########################################################################
