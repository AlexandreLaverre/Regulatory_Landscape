#########################################################################################################################

source("parameters.R")

#########################################################################################################################

for(ref_sp in c("human", "mouse")){
  target_sp = setdiff(c("human", "mouse"), ref_sp)

  print(refs_sp)
  
  enhancers = enhancer.datasets[[ref_sp]]
  
  selenh="ENCODE"
  
  load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", ref_sp, ".RData", sep=""))
  
  ## enhancer statistics
  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
  
  enh.stats.obs=enhancer.statistics[[ref_sp]][[selenh]][["original"]]
  enh.stats.sim=enhancer.statistics[[ref_sp]][[selenh]][["simulated"]]
  
  ## enhancer alignment
  align_enhancer_obs=list_align_enh[[selenh]][["enh_align_obs"]]
  align_enhancer_sim=list_align_enh[[selenh]][["enh_align_simul"]]
  
  ## fragment statistics
  load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))
  
  frag.stats.obs=fragment.statistics[[ref_sp]][["original"]]
  frag.stats.sim=fragment.statistics[[ref_sp]][["simulated"]]
  
  ## add nb genes 
  
  frag_align_obs$nb_genes_500kb=frag.stats.obs[frag_align_obs$ID, "nb_genes_500kb"]
  frag_align_simul$nb_genes_500kb=frag.stats.sim[frag_align_simul$ID, "nb_genes_500kb"]
  
  align_enhancer_obs$nb_genes_500kb=enh.stats.obs[align_enhancer_obs$ID, "nb_genes_500kb"]
  align_enhancer_sim$nb_genes_500kb=enh.stats.sim[align_enhancer_sim$ID, "nb_genes_500kb"]
  
  ## add nb genes class
 
  frag_align_obs$class_genes_500kb=cut(frag_align_obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag_align_obs$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))
  
  frag_align_simul$class_genes_500kb=cut(frag_align_simul$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag_align_simul$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))
  
  align_enhancer_obs$class_genes_500kb=cut(align_enhancer_obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(align_enhancer_obs$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))
  
  align_enhancer_sim$class_genes_500kb=cut(align_enhancer_sim$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(align_enhancer_sim$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))

  ## add repeats
  
  frag.stats.obs$pcrepeat=100*frag.stats.obs$repeat_bp/frag.stats.obs$length
  frag.stats.sim$pcrepeat=100*frag.stats.sim$repeat_bp/frag.stats.sim$length

  frag_align_obs$pcrepeat=frag.stats.obs[frag_align_obs$ID, "pcrepeat"]
  frag_align_simul$pcrepeat=frag.stats.sim[frag_align_simul$ID, "pcrepeat"]
  
  enh.stats.obs$pcrepeat=100*enh.stats.obs$repeat_bp/enh.stats.obs$length
  enh.stats.sim$pcrepeat=100*enh.stats.sim$repeat_bp/enh.stats.sim$length
  
  align_enhancer_obs$pcrepeat=enh.stats.obs[align_enhancer_obs$ID, "pcrepeat"]
  align_enhancer_sim$pcrepeat=enh.stats.sim[align_enhancer_sim$ID, "pcrepeat"]

  ## test

  print("tests for fragments")
  splist=colnames(frag_align_obs)[-c(1, ncol(frag_align_obs):(ncol(frag_align_obs)-4))]
  for(sp in splist){
    print(sp)
    print(paste("observed median fraction of contacted fragment aligned", round(median(frag_align_obs[,sp], na.rm=T)*100,2), "%"))
    print(paste("simulated median fraction of contacted fragment aligned", round(median(frag_align_simul[,sp], na.rm=T)*100,2), "%"))
    print(wilcox.test(frag_align_obs[,sp], frag_align_simul[,sp]))
  }

  print("tests for enhancers")

  for(sp in splist){
    print(sp)
    print(paste("observed median fraction of contacted ENCODE aligned", round(median(align_enhancer_obs[,sp], na.rm=T)*100,2), "%"))
    print(paste("simulated median fraction of contacted ENCODE aligned", round(median(align_enhancer_sim[,sp], na.rm=T)*100,2), "%"))
    print(wilcox.test(align_enhancer_obs[,sp], align_enhancer_sim[,sp]))
  }

 for(type in c("restriction fragments", "enhancers")){
  
   if(type=="restriction fragments"){
     data.obs=frag_align_obs
     data.sim=frag_align_simul
   }

   if(type=="enhancers"){
     data.obs=align_enhancer_obs
     data.sim=align_enhancer_sim
   }
   
   print(paste("computing confidence intervals, distance, ", type))
   
   print("observed")
   BC.obs=tapply(100*data.obs[, other_sp], data.obs$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
   mean.val.obs=unlist(lapply(BC.obs, function(x) x[3]))
   ci.low.obs=unlist(lapply(BC.obs, function(x) x[4]))
   ci.high.obs=unlist(lapply(BC.obs, function(x) x[5]))
   
   print("simulated")
   BC.sim=tapply(100*data.sim[, other_sp], data.sim$dist_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
   mean.val.sim=unlist(lapply(BC.sim, function(x) x[3]))
   ci.low.sim=unlist(lapply(BC.sim, function(x) x[4]))
   ci.high.sim=unlist(lapply(BC.sim, function(x) x[5]))
   
   save(list=c("BC.obs", "BC.sim", "mean.val.obs", "mean.val.sim", "ci.low.obs", "ci.high.obs", "ci.low.sim", "ci.high.sim"), file=paste(pathFigures, "RData/data.bootstrap.conservation.distance.",type,".",ref_sp,".RData"))
   
   
#######################################################################################################
   ## gene density
   
   ## all 
   
   print(paste("computing confidence intervals, gene density, ", type))
   
   print("observed")
   BC.obs=tapply(100*data.obs[, other_sp], data.obs$class_genes_500kb, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
   mean.val.obs=unlist(lapply(BC.obs, function(x) x[3]))
   ci.low.obs=unlist(lapply(BC.obs, function(x) x[4]))
  ci.high.obs=unlist(lapply(BC.obs, function(x) x[5]))
   
  print("simulated")
   BC.sim=tapply(100*data.sim[, other_sp], data.sim$class_genes_500kb, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
   mean.val.sim=unlist(lapply(BC.sim, function(x) x[3]))
   ci.low.sim=unlist(lapply(BC.sim, function(x) x[4]))
   ci.high.sim=unlist(lapply(BC.sim, function(x) x[5]))
   
   ## no repeats
   
   print("observed no repeats")
   
   BC.obs.norep=tapply(100*data.obs[which(data.obs$pcrepeat==0), other_sp], data.obs$class_genes_500kb[which(data.obs$pcrepeat==0)], function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
   mean.val.obs.norep=unlist(lapply(BC.obs.norep, function(x) x[3]))
   ci.low.obs.norep=unlist(lapply(BC.obs.norep, function(x) x[4]))
  ci.high.obs.norep=unlist(lapply(BC.obs.norep, function(x) x[5]))
   
   print("simulated no repeats")
   
   BC.sim.norep=tapply(100*data.sim[which(data.sim$pcrepeat==0), other_sp], data.sim$class_genes_500kb[which(data.sim$pcrepeat==0)], function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
   mean.val.sim.norep=unlist(lapply(BC.sim.norep, function(x) x[3]))
   ci.low.sim.norep=unlist(lapply(BC.sim.norep, function(x) x[4]))
   ci.high.sim.norep=unlist(lapply(BC.sim.norep, function(x) x[5]))
   
   save(list=c("BC.obs", "BC.sim", "mean.val.obs", "mean.val.sim", "ci.low.obs", "ci.high.obs", "ci.low.sim", "ci.high.sim", "BC.obs.norep", "BC.sim.norep", "mean.val.obs.norep", "mean.val.sim.norep", "ci.low.obs.norep", "ci.high.obs.norep", "ci.low.sim.norep", "ci.high.sim.norep"), file=paste(pathFigures, "RData/data.bootstrap.conservation.gene.density.",type,".",ref_sp,".RData"))
   
   print("done")
 }
}
#######################################################################################################
