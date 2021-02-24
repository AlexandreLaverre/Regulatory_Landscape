##########################################################################

objects=ls()

if(!"pathFigures"%in%objects){
  
  source("../main_figures/parameters.R")
  library(bootBCa, lib=pathRlibs)
  
  set.seed(19)
  
  load=T
  prepare=T
}

##########################################################################

if(load){
 ref_sp = "human"

 target_species=c("macaque", "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken") ## other species already done
 
 load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", ref_sp, ".RData", sep=""))

 load=F
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##########################################################################

pdf(paste(pathFigures, "SupplementaryMaterialFigure17.pdf", sep=""), width=6.85, height=10)

## layout

m=matrix(1:8, nrow=4)
layout(m)

###########################################################################

labels=letters[1:length(target_species)]
names(labels)=target_species

par(mar=c(4.1, 4.5, 2.1, 1.5))

for(other_sp in target_species){
  
  nbclasses=length(levels(frag_align_obs$dist_class))
  xpos=1:nbclasses
  
  xlim=c(-0.5, max(xpos)+1)
  
  ## axis position
  class_leg <- c("0", "0.5", "1", "1.5", "2")
  xax=seq(from=0, to=max(xpos)+1, by=10)
  
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

  print("done")
  
  ylim=range(c(ci.low.obs, ci.high.obs, ci.low.simul, ci.high.simul))
  
  dy=diff(ylim)/20
  ylim=ylim+c(-dy, dy)
  
  plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")
  
  ##lines(xpos, mean.val.obs, col=dataset.colors["Original"])
  points(xpos, mean.val.obs, col=dataset.colors["Original"], pch=20)
  segments(xpos, ci.low.obs, xpos, ci.high.obs, col=dataset.colors["Original"])
  
  ##lines(xpos, mean.val.simul, col=dataset.colors["Simulated"])
  points(xpos, mean.val.simul, col=dataset.colors["Simulated"], pch=20)
  segments(xpos, ci.low.simul, xpos, ci.high.simul, col=dataset.colors["Simulated"])
  
  axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)
  mtext("distance to baits (Mb)", side=1, line=2.2, cex=0.8)
  
  axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
  mtext("% aligned sequence", side=2, line=3, cex=0.8)
  
  mtext(paste(ref_sp, " vs. ", other_sp,sep=""), side=3, cex=0.8)
  
  mtext(labels[other_sp], side=3, line=1, at=-7.5, font=2, cex=1.2)
}

###########################################################################


legend("topleft", legend = c("PCHi-C data", "simulated data"), col=dataset.colors, lty=1, bty='n',  inset=c(0, 0.01), xpd=NA, cex=1.25)


###########################################################################

dev.off()

###########################################################################
