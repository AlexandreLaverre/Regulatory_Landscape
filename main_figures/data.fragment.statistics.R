###########################################################################

library(data.table)

source("parameters.R")

pathStats=paste(pathFinalData, "SupplementaryDataset5/", sep="")
pathFragments=paste(pathFinalData, "SupplementaryDataset1/", sep="")

###########################################################################

fragment.statistics=list()

for(sp in c("human", "mouse")){

  print(sp)
  
  obs=fread(paste(pathStats, sp, "/statistics_contacted_sequence_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  sim=fread(paste(pathStats, sp, "/statistics_contacted_sequence_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  
  class(obs)<-"data.frame"
  class(sim)<-"data.frame"

  ## create fragment identifier
  
  obs$ID <-  do.call(paste,c(obs[c("chr","start","end")],sep=":"))
  sim$ID <-  do.call(paste,c(sim[c("chr","start","end")],sep=":"))
  
  ## rownames
  rownames(obs) <- obs$ID
  rownames(sim) <- sim$ID

  ## total covered length

  obs$length=obs$end-obs$start+1
  sim$length=sim$end-sim$start+1

  print(paste("total covered length, observed, before filtering: ",sum(obs$length)))
  print(paste("total covered length, simulated, before filtering: ",sum(sim$length)))
  
  ## remove outlier fragments

  outlier.frag <- read.table(paste(pathFragments, sp, "/aberrant_fragments.txt", sep=""), h=T, stringsAsFactors=F)
 
  print(paste(length(which(outlier.frag$fragment%in%obs$ID)), "outlier fragments in observed data"))
  print(paste(length(which(outlier.frag$fragment%in%sim$ID)), "outlier fragments in simulated data"))
  
  obs <- obs[which(!obs$ID%in%outlier.frag$fragment),]
  sim <- sim[which(!sim$ID%in%outlier.frag$fragment),]
                   
  ## unbaited only
  
  obs <- obs[which(obs$baited == "unbaited"),]
  sim <- sim[which(sim$baited == "unbaited"),]

  ## size selection
 
  obs <- obs[which(obs$length>=minFragmentSize & obs$length<=maxFragmentSize),]
  sim <- sim[which(sim$length>=minFragmentSize & sim$length<=maxFragmentSize),]
  
  ## select fragments depending on number of BLAT hits
   
  obs <- obs[which(obs$BLAT_match > minBLAT & obs$BLAT_match < maxBLAT),] 
  sim <- sim[which(sim$BLAT_match > minBLAT & sim$BLAT_match < maxBLAT),]
  
  ## remove fragment with low read coverage
  coverage = fread(paste(pathFinalData, "SupplementaryDataset1", sp, "reads.coverage.txt", sep="/"), h=T, select=c("ID", "SumReads", "MaxReads"))

  treshold = quantile(coverage$MaxReads, probs=0.1) 
  covered.frag = coverage[which(coverage$MaxReads>treshold),]$ID 
  
  print(paste("coverage threshold: ",treshold, " reads minimum in a sample"))
  
  obs <- obs[which(obs$ID %in% covered.frag),]
  sim <- sim[which(sim$ID %in% covered.frag),]
  
  ## save results

  print(paste("total covered length, observed, after filtering: ",sum(obs$length)))
  print(paste("total covered length, simulated, after filtering: ",sum(sim$length)))
  
  fragment.statistics[[sp]]=list("original"=obs, "simulated"=sim)
}

###########################################################################

save(list=c("fragment.statistics"), file=paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))

###########################################################################
