#############################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

#############################################################################################
library(plyr)

if(load){
  load(paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load=F
}

#############################################################################################

distance.bait.obs.all = list()
distance.bait.sim.all = list()
distance.bait.obs.sim = list()

mean.dist.bait.obs = list()
mean.dist.bait.sim = list()

sd.dist.bait.obs = list()
sd.dist.bait.sim  = list()

distrib.all = list()
example = list()

for (sp in c("human", "mouse")){
  print(sp)
  all.obs = observed.contacts[[sp]]
  all.sim = simulated.contacts[[sp]]
  
  all.obs$dist.bin = cut(all.obs$distance, breaks=seq(from=min(all.obs$distance), to=max(all.obs$distance), by=5000), include.lowest=T)
  all.sim$dist.bin = cut(all.sim$distance, breaks=seq(from=min(all.sim$distance), to=max(all.sim$distance), by=5000), include.lowest=T)
  
  for(sample in sampleinfo[[sp]]$Sample.ID){
    print(sample)
    
    obs = all.obs[which(!is.na(all.obs[sample])),]
    sim = all.sim[which(!is.na(all.sim[sample])),]
    
    ######################################################################################
    ## Density of all interactions by distance bin
    
    dist.obs = t(as.matrix(table(obs$dist.bin)))
    dist.sim = t(as.matrix(table(sim$dist.bin)))
    
    distrib.obs = dist.obs/sum(dist.obs)
    distrib.sim = dist.sim/sum(dist.sim)
    
    distrib.all[[sp]]$density = distrib.obs[1,]
    distrib.all[[sp]]$N = sum(dist.obs)
    
    ######################################################################################
    ## Density of contact by bait and by distance bin 
    
    dist.by.bait.obs = as.data.frame(do.call(rbind, tapply(obs$dist.bin, factor(obs$id_bait), table)))
    dist.by.bait.sim = as.data.frame(do.call(rbind, tapply(sim$dist.bin, factor(sim$id_bait), table)))
    
    # filter on minimum contact number
    dist.by.bait.obs = dist.by.bait.obs[which(apply(dist.by.bait.obs,1,sum) > 5),]
    dist.by.bait.sim = dist.by.bait.sim[which(apply(dist.by.bait.sim,1,sum) > 5),]
    
    common = intersect(rownames(dist.by.bait.obs), rownames(dist.by.bait.sim))
    dist.by.bait.obs = dist.by.bait.obs[common,]
    dist.by.bait.sim = dist.by.bait.sim[common,]
    
    # frequency to probability
    dist.by.bait.obs.prob=dist.by.bait.obs/apply(dist.by.bait.obs,1,sum)
    dist.by.bait.sim.prob=dist.by.bait.sim/apply(dist.by.bait.sim,1,sum)
    
    ######################################################################################
    #### Calculate Euclidean distance
    dist.by.bait.obs.prob=as.matrix(dist.by.bait.obs.prob)
    dist.by.bait.sim.prob=as.matrix(dist.by.bait.sim.prob)
    
    # Bait obs vs all obs
    distance.bait.obs.all[[sp]][[sample]] <- unlist(lapply(1:dim(dist.by.bait.obs.prob)[1], function(x) {
      y=as.numeric(dist.by.bait.obs.prob[x,]); z=as.numeric(distrib.obs[1,]);
      return(sqrt(sum((y-z)^2)))}))
    
    names(distance.bait.obs.all[[sp]][[sample]])=rownames(dist.by.bait.obs.prob)
    
    # Bait sim vs all sim
    distance.bait.sim.all[[sp]][[sample]] <- unlist(lapply(1:dim(dist.by.bait.sim.prob)[1], function(x) {
      y=as.numeric(dist.by.bait.sim.prob[x,]); z=as.numeric(distrib.obs[1,]);
      return(sqrt(sum((y-z)^2)))}))
    
    names(distance.bait.sim.all[[sp]][[sample]])=rownames(dist.by.bait.sim.prob)
    
    # Bait obs vs bait simul
    distance.bait.obs.sim[[sp]][[sample]] <- unlist(lapply(1:dim(dist.by.bait.obs.prob)[1], function(x) {
      y=as.numeric(dist.by.bait.obs.prob[x,]); z=as.numeric(dist.by.bait.sim.prob[x,]);
      return(sqrt(sum((y-z)^2)))}))

    names(distance.bait.obs.sim[[sp]][[sample]])=rownames(dist.by.bait.obs.prob)
    
  }
  
  ######################################################################################
  # Get median of all samples
  
  mean.dist.bait.obs[[sp]]=unlist(lapply(sampleinfo[[sp]]$Sample.ID, function(x) mean(distance.bait.obs.all[[sp]][[x]])))
  mean.dist.bait.sim[[sp]]=unlist(lapply(sampleinfo[[sp]]$Sample.ID, function(x) mean(distance.bait.sim.all[[sp]][[x]])))

  sd.dist.bait.obs[[sp]]=unlist(lapply(sampleinfo[[sp]]$Sample.ID, function(x) sd(distance.bait.obs.all[[sp]][[x]])))
  sd.dist.bait.sim[[sp]]=unlist(lapply(sampleinfo[[sp]]$Sample.ID, function(x) sd(distance.bait.sim.all[[sp]][[x]])))

  ######################################################################################
  # Density of example
  if (sp == "human"){bait = "chr10:95692658:95697241"}else{bait="chr11:71846630:71849155"}
  
  example[[sp]]$distrib.obs.=as.numeric(dist.by.bait.obs.prob[bait,])
  example[[sp]]$distrib.sim=as.numeric(dist.by.bait.sim.prob[bait,])
  
  example[[sp]]$Nobs=sum(dist.by.bait.obs[bait,])
  example[[sp]]$Nsim=sum(dist.by.bait.sim[bait,])
  example[[sp]]$dist.obs=distance.bait.obs.all[[sp]][[sample]][bait]
  example[[sp]]$dist.sim=distance.bait.sim.all[[sp]][[sample]][bait]
}

######################################################################################

save(list=c("distance.bait.obs.all", "distance.bait.sim.all", "distance.bait.obs.sim",
            "mean.dist.bait.obs", "mean.dist.bait.sim", 
            "sd.dist.bait.obs", "sd.dist.bait.sim",
            "example", "distrib.all"), file=paste(pathFigures, "RData/data.bait.euclidean.distances.RData",sep=""))

######################################################################################
