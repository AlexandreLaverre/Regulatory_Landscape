#######################################################################################

source("parameters.R")

pathTables=paste(pathFinalData, "SupplementaryTables/", sep="")

#######################################################################################

suptable1=read.table(paste(pathTables, "SupplementaryTable1.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

#######################################################################################

sampleinfo=list()

for(sp in c("human", "mouse")){
  
  this.info=suptable1[which(suptable1$Species==sp),]

  sampleinfo[[sp]]=this.info
}

#######################################################################################

save(list=c("sampleinfo"), file=paste(pathFigures, "RData/data.sample.info.RData",sep=""))

#######################################################################################
