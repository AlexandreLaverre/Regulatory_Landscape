#################################################################################
library(plyr)

source("parameters.R")

load(paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))

#################################################################################

sp="human"
obs=observed.contacts[[sp]]
sim=simulated.contacts[[sp]]

obs = obs[which(!is.na(obs$Mac1)),]
sim = sim[which(!is.na(sim$Mac1)),]

obs$dist.bin = cut(obs$distance, breaks=seq(from=min(obs$distance), to=max(obs$distance), by=5000), include.lowest=T)
sim$dist.bin = cut(sim$distance, breaks=seq(from=min(sim$distance), to=max(sim$distance), by=5000), include.lowest=T)

######################################################################################
## All interactions 

distrib.obs = t(as.matrix(table(obs$dist.bin)))
distrib.sim = t(as.matrix(table(sim$dist.bin)))

distrib.obs=distrib.obs/sum(distrib.obs)
distrib.sim=distrib.sim/sum(distrib.sim)

par(mfrow=c(2,1))
par(mai=c(1,0.8,0.5,0.5))
plot(distrib.obs[1,], type="l", xlab="Distance bin", ylab="Probability", main="Distribution of all interactions in Mac1")
lines(distrib.sim[1,], col="red")

legend("topright", col=c("black", "red"), legend=c("obs", "sim"), bty="n", lty=1)

eucl.dist = signif(sqrt(sum((distrib.obs[1,]-distrib.sim[1,])^2)), digits=2)
mtext(paste0("Euclidean distance = ", eucl.dist), line=-2)

######################################################################################
## Distribution of contact by distance bin and by bait

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

# Euclidean distance
dist.by.bait.obs.prob=as.matrix(dist.by.bait.obs.prob)
dist.by.bait.sim.prob=as.matrix(dist.by.bait.sim.prob)

distance <- unlist(lapply(1:dim(dist.by.bait.obs.prob)[1], function(x) {
  y=as.numeric(dist.by.bait.obs.prob[x,]); z=as.numeric(dist.by.bait.sim.prob[x,]);
  return(sqrt(sum((y-z)^2)))}))

names(distance)=rownames(dist.by.bait.obs.prob)

hist(distance, breaks=100, xlab="Euclidean Distance between observed and simulated", main=paste0("N = ", length(distance), " baits with at least 5 contacts"))

######################################################################################
## Plot some example
bait_med = "chr1:212351400:212358598" # bait with median number of contact
bait_large = "chr1:209826891:209830217" # bait with large number of contact
bait_dist = "chrX:47076465:47083107" # bait with maximum euclidean distance obs-sim
baits = c(bait_med, bait_large, bait_dist)

for (bait in baits){
  distrib.obs=as.numeric(dist.by.bait.obs.prob[bait,])
  distrib.sim=as.numeric(dist.by.bait.sim.prob[bait,])
  
  Nobs=sum(dist.by.bait.obs[bait,])
  Nsim=sum(dist.by.bait.sim[bait,])
  
  plot(distrib.obs, type="l", xlab="Distance bin", ylab="Probability", main="Bait distribution")
  lines(distrib.sim, col="red")
  
  legend("topright", col=c("black", "red"), legend=c(paste0("obs N=",Nobs),paste0("sim N=",Nsim)), bty="n", lty=1)
  
  mtext(paste0("Euclidean distance = ", signif(distance[bait], digits=2)))

}
######################################################################################




