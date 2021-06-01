
source("../main_figures/parameters.R")

set.seed(19)

library(bootBCa, lib=pathRlibs)

## install.packages("bootBCa", repos="http://R-Forge.R-project.org")

#########################################################################

for(ref in c("human", "mouse")){
  tg=setdiff(c("human", "mouse"), ref)
  
  enh="ENCODE"

  ######################################################################

  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))

  enh.stats.obs=enhancer.statistics[[ref]][[enh]][["original"]]
  enh.stats.sim=enhancer.statistics[[ref]][[enh]][["simulated"]]

  load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref,"2", tg,".RData", sep=""))

  enh.stats.obs$pcungapped=100*pcungapped[rownames(enh.stats.obs)]
  enh.stats.sim$pcungapped=100*pcungapped[rownames(enh.stats.sim)]

  load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))

  frag.stats.obs=fragment.statistics[[ref]][["original"]]
  frag.stats.sim=fragment.statistics[[ref]][["simulated"]]
  
  load(paste(pathFigures, "RData/data.sequence.conservation.fragments.",ref,"2", tg,".RData", sep=""))

  frag.stats.obs$pcungapped=100*pcungapped[rownames(frag.stats.obs)]
  frag.stats.sim$pcungapped=100*pcungapped[rownames(frag.stats.sim)]
  
   ##########################################################################

  ## fragments

  class.frag.nbgenes.obs=cut(frag.stats.obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag.stats.obs$nb_genes_500kb)), include.lowest=T, labels=c("[0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]", ">30"))
  class.frag.nbgenes.sim=cut(frag.stats.sim$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag.stats.sim$nb_genes_500kb)), include.lowest=T, labels=c("[0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]", ">30"))

  ## pc repeats

  frag.stats.obs$pcrepeat=100*frag.stats.obs$repeat_bp/frag.stats.obs$length
  frag.stats.sim$pcrepeat=100*frag.stats.sim$repeat_bp/frag.stats.sim$length
  
  print("computing confidence intervals, repeats vs nb genes, fragments")
  
  BC.obs=tapply(frag.stats.obs$pcrepeat,  class.frag.nbgenes.obs, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.rep.frag.obs=unlist(lapply(BC.obs, function(x) x[3]))
  names(mean.rep.frag.obs)=levels(class.frag.nbgenes.obs)
  
  ci.rep.frag.low.obs=unlist(lapply(BC.obs, function(x) x[4]))
  ci.rep.frag.high.obs=unlist(lapply(BC.obs, function(x) x[5]))

  BC.sim=tapply(frag.stats.sim$pcrepeat,  class.frag.nbgenes.sim, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.rep.frag.sim=unlist(lapply(BC.sim, function(x) x[3]))
  names(mean.rep.frag.sim)=levels(class.frag.nbgenes.sim)
  
  ci.rep.frag.low.sim=unlist(lapply(BC.sim, function(x) x[4]))
  ci.rep.frag.high.sim=unlist(lapply(BC.sim, function(x) x[5]))

  print("done")

  ## conservation as a function of the proportion of repeats

  frag.stats.obs$repeat_class=cut(frag.stats.obs$pcrepeat, breaks=seq(from=0, to=100, length=6), include.lowest=T)
  frag.stats.sim$repeat_class=cut(frag.stats.sim$pcrepeat, breaks=seq(from=0, to=100, length=6), include.lowest=T)

  print("computing confidence intervals, conservation vs. pc repeats, fragments")

  BC.obs=tapply(frag.stats.obs$pcungapped, frag.stats.obs$repeat_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.cons.repclass.frag.obs=unlist(lapply(BC.obs, function(x) x[3]))
  names(mean.cons.repclass.frag.obs)=levels(frag.stats.obs$repeat_class)
  
  ci.low.cons.repclass.frag.obs=unlist(lapply(BC.obs, function(x) x[4]))
  ci.high.cons.repclass.frag.obs=unlist(lapply(BC.obs, function(x) x[5]))
  
  BC.sim=tapply(frag.stats.sim$pcungapped, frag.stats.sim$repeat_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.cons.repclass.frag.sim=unlist(lapply(BC.sim, function(x) x[3]))
  names(mean.cons.repclass.frag.sim)=levels(frag.stats.sim$repeat_class)
  
  ci.low.cons.repclass.frag.sim=unlist(lapply(BC.sim, function(x) x[4]))
  ci.high.cons.repclass.frag.sim=unlist(lapply(BC.sim, function(x) x[5]))
  print("done")
  
  ## enhancers
  
  class.enh.nbgenes.obs=cut(enh.stats.obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(enh.stats.obs$nb_genes_500kb)), include.lowest=T, labels=c("[0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]", ">30"))
  class.enh.nbgenes.sim=cut(enh.stats.sim$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(enh.stats.sim$nb_genes_500kb)), include.lowest=T, labels=c("[0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]", ">30"))

 
  ## pc repeats

  enh.stats.obs$pcrepeat=100*enh.stats.obs$repeat_bp/enh.stats.obs$length
  enh.stats.sim$pcrepeat=100*enh.stats.sim$repeat_bp/enh.stats.sim$length
  
  print("computing confidence intervals, pc repeats vs. genes, enhancers")
  
  BC.obs=tapply(enh.stats.obs$pcrepeat,  class.enh.nbgenes.obs, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.rep.enh.obs=unlist(lapply(BC.obs, function(x) x[3]))
  names(mean.rep.enh.obs)=levels(class.enh.nbgenes.obs)
        
  ci.rep.enh.low.obs=unlist(lapply(BC.obs, function(x) x[4]))
  ci.rep.enh.high.obs=unlist(lapply(BC.obs, function(x) x[5]))
  
  BC.sim=tapply(enh.stats.sim$pcrepeat,  class.enh.nbgenes.sim, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.rep.enh.sim=unlist(lapply(BC.sim, function(x) x[3]))
  names(mean.rep.enh.sim)=levels(class.enh.nbgenes.sim)
  
  ci.rep.enh.low.sim=unlist(lapply(BC.sim, function(x) x[4]))
  ci.rep.enh.high.sim=unlist(lapply(BC.sim, function(x) x[5]))

  print("done")

  ## conservation as a function of the proportion of repeats
  
  enh.stats.obs$repeat_class=cut(enh.stats.obs$pcrepeat, breaks=seq(from=0, to=100, length=6), include.lowest=T)
  enh.stats.sim$repeat_class=cut(enh.stats.sim$pcrepeat, breaks=seq(from=0, to=100, length=6), include.lowest=T)

  print("computing confidence intervals, conservation vs. pc repeats, fragments")
  
  BC.obs=tapply(enh.stats.obs$pcungapped,  enh.stats.obs$repeat_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.cons.repclass.enh.obs=unlist(lapply(BC.obs, function(x) x[3]))
  names(mean.cons.repclass.enh.obs)=levels(enh.stats.obs$repeat_class)
  
  ci.low.cons.repclass.enh.obs=unlist(lapply(BC.obs, function(x) x[4]))
  ci.high.cons.repclass.enh.obs=unlist(lapply(BC.obs, function(x) x[5]))
  
  BC.sim=tapply(enh.stats.sim$pcungapped,  enh.stats.sim$repeat_class, function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
  mean.cons.repclass.enh.sim=unlist(lapply(BC.sim, function(x) x[3]))
  names(mean.cons.repclass.enh.sim)=levels(enh.stats.sim$repeat_class)
        
  ci.low.cons.repclass.enh.sim=unlist(lapply(BC.sim, function(x) x[4]))
  ci.high.cons.repclass.enh.sim=unlist(lapply(BC.sim, function(x) x[5]))

  print("done")

  ## save objects

  save(list=ls(), file=paste(pathFigures, "RData/data.bootstrap.sequence.features.",ref,"RData",sep=""))
}

##########################################################################
