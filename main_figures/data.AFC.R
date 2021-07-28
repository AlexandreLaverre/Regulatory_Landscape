########################################################################

source("parameters.R")

########################################################################

load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

########################################################################

library(ade4)

subsample = TRUE

if (subsample){sample.suffix=".subsample"}else{sample.suffix=""}
  
########################################################################

data.AFC=list()

for(sp in c("human", "mouse")){
  info=sampleinfo[[sp]]
  samples=info$Sample.ID

  data_AFC <- observed.contacts[[sp]]
  data_AFC <- data_AFC[, samples] 
  
  # Sub-sampling to minimum number of interactions
  if (subsample){
    min.nb.interaction = min(apply(data_AFC, 2, function(x) nrow(data_AFC[which(!is.na(x)),])))
    
    for (sample in samples){
      treshold.value =  min(head(sort(data_AFC[[sample]] ,decreasing=TRUE), n=min.nb.interaction))
      data_AFC[which(data_AFC[[sample]] < treshold.value), sample] = NA
    }
  }

  data_AFC [!is.na(data_AFC)] <- 1
  data_AFC [is.na(data_AFC)] <- 0
  data_AFC <- data.frame(t(data_AFC)) # row = samples

  this.AFC <- dudi.coa(data_AFC, scannf=F, nf=3)

  dist.AFC=dist.dudi(this.AFC) 
  hclust.AFC=hclust(dist.AFC, "ward")
  
  data.AFC[[sp]]=list("AFC"=this.AFC, "hclust.AFC"=hclust.AFC)
}

########################################################################

save(data.AFC, file=paste(pathFigures, "RData/data.AFC", sample.suffix, ".RData", sep=""))

########################################################################
