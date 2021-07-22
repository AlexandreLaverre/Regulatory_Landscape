##########################################################################

source("../main_figures/parameters.R")
library(bootBCa, lib=pathRlibs)

set.seed(19)

for(ref_sp in c( "human", "mouse")){
  
  print(ref_sp)
  load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.and.phyloP.", ref_sp, ".RData", sep=""))
  
  ### enhancer
  #for(selenh in c(enhancer.datasets[[ref_sp]])){
  selenh="ENCODE"
  align_enhancer_obs=list_align_enh[[selenh]][["enh_align_obs"]]
  align_enhancer_simul=list_align_enh[[selenh]][["enh_align_simul"]]
  
  print("computing confidence intervals for ENCODE")
  
  print("observed")
  BC.obs=tapply(align_enhancer_obs$phyloP_score, align_enhancer_obs$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.val.obs=unlist(lapply(BC.obs, function(x) x[3]))
  ci.low.obs=unlist(lapply(BC.obs, function(x) x[4]))
  ci.high.obs=unlist(lapply(BC.obs, function(x) x[5]))
  
  print("simulated")
  BC.sim=tapply(align_enhancer_simul$phyloP_score, align_enhancer_simul$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.val.simul=unlist(lapply(BC.sim, function(x) x[3]))
  ci.low.simul=unlist(lapply(BC.sim, function(x) x[4]))
  ci.high.simul=unlist(lapply(BC.sim, function(x) x[5]))
  
  save(list=c("BC.obs", "mean.val.obs", "ci.low.obs", "ci.high.obs", "BC.sim", "mean.val.simul", "ci.low.simul", "ci.high.simul"), 
       file=paste(pathFigures, "RData/data.bootstrap.",selenh,".conservation.distance.",ref_sp,".phyloP.RData",sep=""))
  #}
  
  ### restriction fragments
  print("computing confidence intervals for restriction fragments")
  
  print("observed")
  BC.obs=tapply(frag_align_obs$phyloP_score, frag_align_obs$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.val.obs=unlist(lapply(BC.obs, function(x) x[3]))
  ci.low.obs=unlist(lapply(BC.obs, function(x) x[4]))
  ci.high.obs=unlist(lapply(BC.obs, function(x) x[5]))
  
  print("simulated")
  BC.sim=tapply(frag_align_simul$phyloP_score, frag_align_simul$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.val.simul=unlist(lapply(BC.sim, function(x) x[3]))
  ci.low.simul=unlist(lapply(BC.sim, function(x) x[4]))
  ci.high.simul=unlist(lapply(BC.sim, function(x) x[5]))
  
  
  save(list=c("BC.obs", "mean.val.obs", "ci.low.obs", "ci.high.obs", "BC.sim", "mean.val.simul", "ci.low.simul", "ci.high.simul"),
       file=paste(pathFigures, "RData/data.bootstrap.fragments.conservation.distance.", ref_sp,".phyloP.RData",sep="")) 
  
}

####################################################################################################3
