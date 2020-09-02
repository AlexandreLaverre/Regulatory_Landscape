#######################################################################################

source("parameters.R")

pathTables=paste(pathFinalData, "SupplementaryTables/", sep="")

#######################################################################################

separators=c(",", "\t")
names(separators)=c("human", "mouse")

sampleinfo=list()

for(sp in c("human", "mouse")){
  
  this.info=read.table(paste(pathTables, sp, "_samples_informations.csv",sep=""), sep=separators[sp], h=T, stringsAsFactors=F)

  sampleinfo[[sp]]=this.info
}

#######################################################################################

save(list=c("sampleinfo"), file=paste(pathFigures, "RData/data.sample.info.RData",sep=""))

#######################################################################################
