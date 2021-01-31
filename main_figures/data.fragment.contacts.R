#################################################################################

library(data.table)

source("parameters.R")

load(paste(pathFigures, "RData/data.fragment.statistics.RData",sep=""))

#################################################################################

observed.contacts=list()
simulated.contacts=list()

for(sp in c("human", "mouse")){

  ## fragment statistics
  
  frag.stats.obs <- fragment.statistics[[sp]][["original"]]
  frag.stats.sim <- fragment.statistics[[sp]][["simulated"]]
  
  ## interactions
  
  obs <- fread(paste(pathFinalData, "SupplementaryDataset1/", sp, "/all_interactions.txt", sep=""), header=T, stringsAsFactors=F, sep="\t")
  sim <- fread(paste(pathFinalData, "SupplementaryDataset2/", sp, "/simulated_all_interactions.txt", sep=""), header=T, stringsAsFactors=F, sep="\t")

  class(obs)<- "data.frame"
  class(sim)<- "data.frame"
  
  ## select interactions in cis
  obs=obs[which(obs$chr_bait==obs$chr),]
  sim=sim[which(sim$chr_bait==sim$chr),]
  
  ## select bait-other interactions
  obs=obs[which(obs$type=="unbaited"),]
  sim=sim[which(sim$type=="unbaited"),]
  
  ## select interactions in the accepted distance range
  obs=obs[which(obs$distance>=minDistance & obs$distance<=maxDistance),]
  sim=sim[which(sim$distance>=minDistance & sim$distance<=maxDistance),]

  ## discard interactions with unusual restriction fragments length
  obs$bait_length = obs$end_bait-obs$start_bait+1
  obs$contacted_length = obs$end-obs$start+1
  sim$bait_length = sim$end_bait-sim$start_bait+1
  sim$contacted_length = sim$end-sim$start+1
  
  obs=obs[which(obs$bait_length>=minFragmentSize & obs$bait_length<=maxFragmentSize & obs$contacted_length>=minFragmentSize & obs$contacted_length<=maxFragmentSize),]
  sim=sim[which(sim$bait_length>=minFragmentSize & sim$bait_length<=maxFragmentSize & sim$contacted_length>=minFragmentSize & sim$contacted_length<=maxFragmentSize),]
    
  ## bait and fragment id
  obs$id_bait=paste(obs$chr_bait, obs$start_bait, obs$end_bait, sep=":")
  obs$id_frag=paste(obs$chr, obs$start, obs$end, sep=":")
  
  sim$id_bait=paste(sim$chr_bait, sim$start_bait, sim$end_bait, sep=":")
  sim$id_frag=paste(sim$chr, sim$start, sim$end, sep=":")
  
  ## keep only previously filtered fragments

  obs=obs[which(obs$id_frag%in%frag.stats.obs$ID),]
  sim=sim[which(sim$id_frag%in%frag.stats.sim$ID),]

  ## save results
  
  observed.contacts[[sp]]=obs
  simulated.contacts[[sp]]=sim
}

#################################################################################

save(list=c("observed.contacts", "simulated.contacts"), file=paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))

#################################################################################



