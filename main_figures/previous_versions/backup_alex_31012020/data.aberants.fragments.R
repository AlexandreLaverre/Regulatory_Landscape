#################################################################################

library(data.table)

source("parameters.R")
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

#################################################################################

for(sp in c("human", "mouse")){
  obs <- fread(paste(pathFinalData, "SupplementaryDataset1/", sp, "/all_interactions.txt", sep=""), header=T, stringsAsFactors=F, sep="\t")

  class(obs)<- "data.frame"

  ## select interactions in cis
  obs=obs[which(obs$chr_bait==obs$chr),]
  obs=obs[which(obs$distance>=minDistance & obs$distance<=maxDistance),]
  
  ## bait and fragment id
  obs$id_bait=paste(obs$chr_bait, obs$start_bait, obs$end_bait, sep=":")
  obs$id_frag=paste(obs$chr, obs$start, obs$end, sep=":")
  
  ## CHICAGO median score
  samples=sampleinfo[[sp]]$Sample.ID
  obs$median_score <- apply(obs[samples], 1, function(x) median(x, na.rm=T))

  ## discard abberant fragments 
  nb_bait <- tapply(obs$id_bait, as.factor(obs$id_frag), function(x) length(x))
  mean_score <- tapply(obs$median_score, as.factor(obs$id_frag), function(x) mean(x))
  med_dist <- tapply(obs$dist, as.factor(obs$id_frag), function(x) median(x))
  
  frag_stats = data.frame("nb.bait"=nb_bait, "mean.CHICAGO.score"=mean_score, "median.distance"=med_dist)
  aberant_frags = frag_stats[which(frag_stats$nb.bait > 15 & frag_stats$mean.CHICAGO.score>10 & frag_stats$median.distance>1000000),]
  aberant_frags = data.frame("fragment"=rownames(aberant_frags), aberant_frags)
  
  #################################################################################
  
  write.table(aberant_frags, file=paste(pathFinalData, "SupplementaryDataset1/", sp, "/aberants.fragments.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
  
  #################################################################################
  
}








