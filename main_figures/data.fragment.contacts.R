#################################################################################

library(data.table)

source("parameters.R")

#################################################################################

observed.contacts=list()
simulated.contacts=list()

for(sp in c("human", "mouse")){
  obs <- fread(paste(pathFinalData, "SupplementaryDataset1/", sp, "/all_interactions.txt", sep=""), header=T, stringsAsFactors=F, sep="\t")
  sim <- fread(paste(pathFinalData, "SupplementaryDataset2/", sp, "/simulated_all_interactions.txt", sep=""), header=T, stringsAsFactors=F, sep="\t")
  
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
  
  obs=obs[which(obs$bait_length>=250 & obs$bait_length<=20000 & obs$contacted_length>=250 & obs$contacted_length<=20000),]
  sim=sim[which(sim$bait_length>=250 & sim$bait_length<=20000 & sim$contacted_length>=250 & sim$contacted_length<=20000),]
    
  ## bait and fragment id

  obs$id_bait=paste(obs$chr_bait, obs$start_bait, obs$end_bait, sep=":")
  obs$id_frag=paste(obs$chr, obs$start, obs$end, sep=":")
  
  sim$id_bait=paste(sim$chr_bait, sim$start_bait, sim$end_bait, sep=":")
  sim$id_frag=paste(sim$chr, sim$start, sim$end, sep=":")
  
  observed.contacts[[sp]]=obs
  simulated.contacts[[sp]]=sim
}

#################################################################################

save(list=c("observed.contacts", "simulated.contacts"), file=paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))

#################################################################################



