##########################################################################

objects=ls()

if(!"pathFigures"%in%objects){

 source("../main_figures/parameters.R")

 load=T
 prepare=T
}

##########################################################################

if(load){
 ref_sp = "mouse"

 target_species=c("rat", "rabbit", "human", "macaque", "dog", "cow", "elephant", "opossum", "chicken") 
 
 load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", ref_sp, ".RData", sep=""))

 selenh="ENCODE"

 align_enhancer_obs=list_align_enh[[selenh]][["enh_align_obs"]]
 align_enhancer_simul=list_align_enh[[selenh]][["enh_align_simul"]]
 
 load=F
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##########################################################################

pdf(paste(pathFigures, "SupplementaryFigure17.pdf", sep=""), width=6.85, height=11)

## layout

m=matrix(1:10, nrow=5)
layout(m)

###########################################################################

labels=letters[1:length(target_species)]
names(labels)=target_species

par(mar=c(4.1, 4.5, 2.1, 1.5))

for(other_sp in target_species){
  
  nbclasses=length(levels(align_enhancer_obs$dist_class))
  xpos=1:nbclasses
  
  xlim=c(-0.5, max(xpos)+1)
  
  ## axis position
  class_leg <- c("0", "0.5", "1", "1.5", "2")
  xax=seq(from=0, to=max(xpos)+1, by=10)
  
  mean.val.obs=tapply(100*align_enhancer_obs[, other_sp], align_enhancer_obs$dist_class, function(x) mean(x, na.rm=T))
  ci.low.obs=tapply(100*align_enhancer_obs[, other_sp], align_enhancer_obs$dist_class, function(x) t.test(x)[["conf.int"]][1])
  ci.high.obs=tapply(100*align_enhancer_obs[, other_sp], align_enhancer_obs$dist_class, function(x) t.test(x)[["conf.int"]][2])
  
  mean.val.simul=tapply(100*align_enhancer_simul[, other_sp], align_enhancer_simul$dist_class, function(x) mean(x, na.rm=T))
  ci.low.simul=tapply(100*align_enhancer_simul[, other_sp], align_enhancer_simul$dist_class, function(x) t.test(x)[["conf.int"]][1])
  ci.high.simul=tapply(100*align_enhancer_simul[, other_sp], align_enhancer_simul$dist_class, function(x) t.test(x)[["conf.int"]][2])
  
  ylim=range(c(ci.low.obs, ci.high.obs, ci.low.simul, ci.high.simul))
  
  dy=diff(ylim)/20
  ylim=ylim+c(-dy, dy)
  
  plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")
  
  lines(xpos, mean.val.obs, col=dataset.colors["Original"])
  segments(xpos, ci.low.obs, xpos, ci.high.obs, col=dataset.colors["Original"])
  
  lines(xpos, mean.val.simul, col=dataset.colors["Simulated"])
  segments(xpos, ci.low.simul, xpos, ci.high.simul, col=dataset.colors["Simulated"])
  
  axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)
  mtext("distance to promoters (Mb)", side=1, line=2.2, cex=0.8)
  
  axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
  mtext("% aligned sequence", side=2, line=3, cex=0.8)
  
  mtext(paste(ref_sp, " vs. ", other_sp,sep=""), side=3, cex=0.8)
  
  mtext(labels[other_sp], side=3, line=0.8, at=-7.5, font=2, cex=1.2)
}

###########################################################################

## empty plot + legend

plot(1, type="n", xlab="", ylab="", axes=F)


legend("topleft", legend = c("PCHi-C data", "simulated data"), col=dataset.colors, lty=1, bty='n',  inset=c(0, 0.01), xpd=NA, cex=1.25)


###########################################################################

dev.off()

###########################################################################
