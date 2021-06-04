##########################################################################

source("../main_figures/parameters.R")
library(bootBCa, lib=pathRlibs)

set.seed(19)

for(ref_sp in c( "human", "mouse"){
  
  if(ref_sp=="human"){
    target_species=c("macaque", "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken") ## other species already done
  }

  if(ref_sp=="mouse"){
    target_species=c("rat", "rabbit", "human", "macaque", "dog", "cow", "elephant", "opossum", "chicken") }
  
  load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", ref_sp, ".RData", sep=""))
  
  for(other_sp in target_species){
    
    print("computing confidence intervals")
    
    print("observed")
    BC.obs=tapply(100*frag_align_obs[, other_sp], frag_align_obs$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    mean.val.obs=unlist(lapply(BC.obs, function(x) x[3]))
    ci.low.obs=unlist(lapply(BC.obs, function(x) x[4]))
    ci.high.obs=unlist(lapply(BC.obs, function(x) x[5]))
    
    print("simulated")
    BC.sim=tapply(100*frag_align_simul[, other_sp], frag_align_simul$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    mean.val.simul=unlist(lapply(BC.sim, function(x) x[3]))
    ci.low.simul=unlist(lapply(BC.sim, function(x) x[4]))
    ci.high.simul=unlist(lapply(BC.sim, function(x) x[5]))
    
    
    save(list=c("BC.obs", "mean.val.obs", "ci.low.obs", "ci.high.obs", "BC.sim", "mean.val.sim", "ci.low.sim", "ci.high.sim"), file=paste(pathFigures, "RData/data.bootstrap.fragments.conservation.distance.",ref_sp,".",other_sp,".RData",sep="")) 
    
  }

}

####################################################################################################3
