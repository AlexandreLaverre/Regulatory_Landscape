#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  
  source("../main_figures/parameters.R")

  ref="human"
  tg=setdiff(c("human", "mouse"), ref)
  
  enh="ENCODE"
}

##########################################################################

if(load){

  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))

  enh.stats.obs=enhancer.statistics[[ref]][[enh]][["original"]]
  enh.stats.sim=enhancer.statistics[[ref]][[enh]][["simulated"]]

  load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref,"2", tg,".RData", sep=""))

  enh.stats.obs$pcungapped=pcungapped[rownames(enh.stats.obs)]
  enh.stats.sim$pcungapped=pcungapped[rownames(enh.stats.sim)]

  load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))

  frag.stats.obs=fragment.statistics[[ref]][["original"]]
  frag.stats.sim=fragment.statistics[[ref]][["simulated"]]
  
  load(paste(pathFigures, "RData/data.sequence.conservation.fragments.",ref,"2", tg,".RData", sep=""))

  frag.stats.obs$pcungapped=pcungapped[rownames(frag.stats.obs)]
  frag.stats.sim$pcungapped=pcungapped[rownames(frag.stats.sim)]
  
  load=FALSE
}

##########################################################################

if(prepare){

  ## fragments

  class.frag.nbgenes.obs=cut(frag.stats.obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag.stats.obs$nb_genes_500kb)), include.lowest=T, labels=c("[0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]", ">30"))
  class.frag.nbgenes.sim=cut(frag.stats.sim$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag.stats.sim$nb_genes_500kb)), include.lowest=T, labels=c("[0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]", ">30"))

  ## conservation
  
  mean.cons.frag.obs=tapply(frag.stats.obs$pcungapped, class.frag.nbgenes.obs, mean, na.rm=T)
  ci.cons.frag.low.obs=tapply(frag.stats.obs$pcungapped, class.frag.nbgenes.obs, function(x) t.test(x)[["conf.int"]][1])
  ci.cons.frag.high.obs=tapply(frag.stats.obs$pcungapped, class.frag.nbgenes.obs, function(x) t.test(x)[["conf.int"]][2])
  
  mean.cons.frag.sim=tapply(frag.stats.sim$pcungapped, class.frag.nbgenes.sim, mean, na.rm=T)
  ci.cons.frag.low.sim=tapply(frag.stats.sim$pcungapped, class.frag.nbgenes.sim, function(x) t.test(x)[["conf.int"]][1])
  ci.cons.frag.high.sim=tapply(frag.stats.sim$pcungapped, class.frag.nbgenes.sim, function(x) t.test(x)[["conf.int"]][2])

  ## pc repeats

  frag.stats.obs$pcrepeat=100*frag.stats.obs$repeat_bp/frag.stats.obs$length
  frag.stats.sim$pcrepeat=100*frag.stats.sim$repeat_bp/frag.stats.sim$length

  mean.rep.frag.obs=tapply(frag.stats.obs$pcrepeat, class.frag.nbgenes.obs, mean, na.rm=T)
  ci.rep.frag.low.obs=tapply(frag.stats.obs$pcrepeat, class.frag.nbgenes.obs, function(x) t.test(x)[["conf.int"]][1])
  ci.rep.frag.high.obs=tapply(frag.stats.obs$pcrepeat, class.frag.nbgenes.obs, function(x) t.test(x)[["conf.int"]][2])

  
  mean.rep.frag.sim=tapply(frag.stats.sim$pcrepeat, class.frag.nbgenes.sim, mean, na.rm=T)
  ci.rep.frag.low.sim=tapply(frag.stats.sim$pcrepeat, class.frag.nbgenes.sim, function(x) t.test(x)[["conf.int"]][1])
  ci.rep.frag.high.sim=tapply(frag.stats.sim$pcrepeat, class.frag.nbgenes.sim, function(x) t.test(x)[["conf.int"]][2])
  
  
  ## enhancers
  
  class.enh.nbgenes.obs=cut(enh.stats.obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(enh.stats.obs$nb_genes_500kb)), include.lowest=T, labels=c("[0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]", ">30"))
  class.enh.nbgenes.sim=cut(enh.stats.sim$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(enh.stats.sim$nb_genes_500kb)), include.lowest=T, labels=c("[0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]", ">30"))

  ## conservation
  
  mean.cons.enh.obs=tapply(enh.stats.obs$pcungapped, class.enh.nbgenes.obs, mean, na.rm=T)
  ci.cons.enh.low.obs=tapply(enh.stats.obs$pcungapped, class.enh.nbgenes.obs, function(x) t.test(x)[["conf.int"]][1])
  ci.cons.enh.high.obs=tapply(enh.stats.obs$pcungapped, class.enh.nbgenes.obs, function(x) t.test(x)[["conf.int"]][2])
  
  mean.cons.enh.sim=tapply(enh.stats.sim$pcungapped, class.enh.nbgenes.sim, mean, na.rm=T)
  ci.cons.enh.low.sim=tapply(enh.stats.sim$pcungapped, class.enh.nbgenes.sim, function(x) t.test(x)[["conf.int"]][1])
  ci.cons.enh.high.sim=tapply(enh.stats.sim$pcungapped, class.enh.nbgenes.sim, function(x) t.test(x)[["conf.int"]][2])

  ## pc repeats

  enh.stats.obs$pcrepeat=100*enh.stats.obs$repeat_bp/enh.stats.obs$length
  enh.stats.sim$pcrepeat=100*enh.stats.sim$repeat_bp/enh.stats.sim$length

  mean.rep.enh.obs=tapply(enh.stats.obs$pcrepeat, class.enh.nbgenes.obs, mean, na.rm=T)
  ci.rep.enh.low.obs=tapply(enh.stats.obs$pcrepeat, class.enh.nbgenes.obs, function(x) t.test(x)[["conf.int"]][1])
  ci.rep.enh.high.obs=tapply(enh.stats.obs$pcrepeat, class.enh.nbgenes.obs, function(x) t.test(x)[["conf.int"]][2])

  
  mean.rep.enh.sim=tapply(enh.stats.sim$pcrepeat, class.enh.nbgenes.sim, mean, na.rm=T)
  ci.rep.enh.low.sim=tapply(enh.stats.sim$pcrepeat, class.enh.nbgenes.sim, function(x) t.test(x)[["conf.int"]][1])
  ci.rep.enh.high.sim=tapply(enh.stats.sim$pcrepeat, class.enh.nbgenes.sim, function(x) t.test(x)[["conf.int"]][2])
  
  prepare=F
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "ExtendedFigure4.pdf", sep=""), width=4.49, height=6)

m=matrix(rep(NA, 2*6), nrow=2)
m[1,]=c(rep(1,2), rep(2, 4))
m[2,]=c(rep(3,2), rep(4, 4))

layout(m)

###########################################################################

labels.genes=c("a", "c")
labels.repeats=c("b", "d")

names(labels.genes)=c("frag", "enh")
names(labels.repeats)=c("frag", "enh")

for(type in c("frag", "enh")){

  ## get objects
  stats.obs=get(paste(type, ".stats.obs", sep=""))
  stats.sim=get(paste(type, ".stats.sim", sep=""))

  mean.cons.obs=get(paste("mean.cons.",type,".obs",sep=""))
  ci.cons.low.obs=get(paste("ci.cons.",type,".low.obs",sep=""))
  ci.cons.high.obs=get(paste("ci.cons.",type,".high.obs",sep=""))
  
  mean.rep.obs=get(paste("mean.rep.",type,".obs",sep=""))
  ci.rep.low.obs=get(paste("ci.rep.",type,".low.obs",sep=""))
  ci.rep.high.obs=get(paste("ci.rep.",type,".high.obs",sep=""))

  mean.cons.sim=get(paste("mean.cons.",type,".sim",sep=""))
  ci.cons.low.sim=get(paste("ci.cons.",type,".low.sim",sep=""))
  ci.cons.high.sim=get(paste("ci.cons.",type,".high.sim",sep=""))
  
  mean.rep.sim=get(paste("mean.rep.",type,".sim",sep=""))
  ci.rep.low.sim=get(paste("ci.rep.",type,".low.sim",sep=""))
  ci.rep.high.sim=get(paste("ci.rep.",type,".high.sim",sep=""))
  
  ###########################################################################

  ## nb genes
  
  par(mar=c(6.75, 3.75, 2.1, 1.1))

  if(ref=="human"){
    ylim=range(c(0, 25))
  }

   if(ref=="mouse"){
    ylim=range(c(0, 40))
  } 

  plot(1, type="n", xlim=c(0.5, 2.5), ylim=ylim, xlab="", ylab="", axes=F)
  
  boxplot(stats.obs$nb_genes_500kb, col="white", border=dataset.colors[["Original"]], add=T, at=1, axes=F, outline=F, boxwex=0.5, notch=T)
  
  boxplot(stats.sim$nb_genes_500kb, pch=20, col="white", border=dataset.colors[["Simulated"]], add=T, at=2, axes=F, outline=F, boxwex=0.5, notch=T)
  
  axis(side=1, mgp=c(3, 0.5,0), at=1:2, labels=rep("", 2))
  axis(side=2, mgp=c(3, 0.75,0))
  
  mtext("number of genes within 500kb", side=2, line=2.5, cex=0.75)
  
  mtext(c("PCHi-C", "simulated"), at=1:2, side=1, line=0.75, cex=0.75, las=2)
  
  mtext(labels.genes[type], side=3, at=-0.55, font=2, line=1)
 
 ###########################################################################

  ## pc repeats

  par(mar=c(6.75, 3.5, 2.1, 1))
  
  xpos=1:length(mean.rep.obs)
  xlim=range(xpos)+c(-0.5, 0.5)

  ylim=range(c(ci.rep.low.obs, ci.rep.high.obs, ci.rep.low.sim, ci.rep.high.sim))
  smally=diff(ylim)/10
  ylim=ylim+c(-smally,smally)
  
  
  smallx=c(-0.1, 0.1)
  
  plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim)
  
  points(xpos+smallx[1], mean.rep.obs, col=dataset.colors["Original"], pch=20)
  segments(xpos+smallx[1], ci.rep.low.obs, xpos+smallx[1], ci.rep.high.obs, col=dataset.colors[["Original"]])
  
  points(xpos+smallx[2], mean.rep.sim, col=dataset.colors["Simulated"], pch=20)
  segments(xpos+smallx[2], ci.rep.low.sim, xpos+smallx[2], ci.rep.high.sim, col=dataset.colors[["Simulated"]])
  
  abline(v=xpos[-length(xpos)]+0.5, lty=3, col="gray40")
  
  axis(side=2, mgp=c(3, 0.75, 0), las=2)
  axis(side=1, mgp=c(3, 0.75, 0), at=xpos, labels=rep("", length(xpos)))
  
  mtext(names(mean.rep.obs), at=xpos, side=1, line=1, las=2, cex=0.75)
  
  mtext("number of genes within 500 kb", side=1, line=4.5, cex=0.75)
  mtext("% repetitive sequence", side=2, line=2.5, cex=0.75)

  mtext(labels.repeats[type], side=3, at=-1, font=2, line=1)

  
  if(type=="frag"){
    legend("topleft", col=dataset.colors, legend = c("PCHi-C data", "simulated data"), box.col="white", bg="white", pch=20, cex=1.1, inset=c(0.01,-0.05), xpd=NA)

    mtext("restriction fragments", side=4, line=-0.25, cex=0.75)
  } else{
    mtext("enhancers", side=4, line=-0.25, cex=0.75)
  }  

}

###########################################################################

dev.off()

###########################################################################
