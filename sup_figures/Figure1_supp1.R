#############################################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

#############################################################################################

if(load){
 
  load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

  load=FALSE
}

#############################################################################################

if(prepare){
  sampleinfo.human=sampleinfo[["human"]]
  sampleinfo.mouse=sampleinfo[["mouse"]]

  exsample.human=sampleinfo.human$Sample.ID[1]
  exsample.mouse=sampleinfo.mouse$Sample.ID[1]

  ## contact data

  obs.human=observed.contacts[["human"]]
  sim.human=simulated.contacts[["human"]]

  obs.mouse=observed.contacts[["mouse"]]
  sim.mouse=simulated.contacts[["mouse"]]

  ## median distance per sample

  median.distance.obs.human=unlist(lapply(sampleinfo.human$Sample.ID, function(x) median(obs.human[which(!is.na(obs.human[,x])),"distance"])))
  median.distance.sim.human=unlist(lapply(sampleinfo.human$Sample.ID, function(x) median(sim.human[which(!is.na(sim.human[,x])),"distance"])))

  median.distance.obs.mouse=unlist(lapply(sampleinfo.mouse$Sample.ID, function(x) median(obs.mouse[which(!is.na(obs.mouse[,x])),"distance"])))
  median.distance.sim.mouse=unlist(lapply(sampleinfo.mouse$Sample.ID, function(x) median(sim.mouse[which(!is.na(sim.mouse[,x])),"distance"])))


  ## mean nb contacts per bait

  mean.nb.contacts.bait.obs.human=unlist(lapply(sampleinfo.human$Sample.ID, function(x) {this.obs=obs.human[which(!is.na(obs.human[,x])),]; return(mean(as.numeric(table(this.obs$id_bait))))}))
  mean.nb.contacts.bait.sim.human=unlist(lapply(sampleinfo.human$Sample.ID, function(x) {this.sim=sim.human[which(!is.na(sim.human[,x])),]; return(mean(as.numeric(table(this.sim$id_bait))))}))
  
  mean.nb.contacts.bait.obs.mouse=unlist(lapply(sampleinfo.mouse$Sample.ID, function(x) {this.obs=obs.mouse[which(!is.na(obs.mouse[,x])),]; return(mean(as.numeric(table(this.obs$id_bait))))}))
  mean.nb.contacts.bait.sim.mouse=unlist(lapply(sampleinfo.mouse$Sample.ID, function(x) {this.sim=sim.mouse[which(!is.na(sim.mouse[,x])),]; return(mean(as.numeric(table(this.sim$id_bait))))}))

  ## mean nb contacts per fragment

  mean.nb.contacts.frag.obs.human=unlist(lapply(sampleinfo.human$Sample.ID, function(x) {this.obs=obs.human[which(!is.na(obs.human[,x])),]; return(mean(as.numeric(table(this.obs$id_frag))))}))
  mean.nb.contacts.frag.sim.human=unlist(lapply(sampleinfo.human$Sample.ID, function(x) {this.sim=sim.human[which(!is.na(sim.human[,x])),]; return(mean(as.numeric(table(this.sim$id_frag))))}))

  mean.nb.contacts.frag.obs.mouse=unlist(lapply(sampleinfo.mouse$Sample.ID, function(x) {this.obs=obs.mouse[which(!is.na(obs.mouse[,x])),]; return(mean(as.numeric(table(this.obs$id_frag))))}))
  mean.nb.contacts.frag.sim.mouse=unlist(lapply(sampleinfo.mouse$Sample.ID, function(x) {this.sim=sim.mouse[which(!is.na(sim.mouse[,x])),]; return(mean(as.numeric(table(this.sim$id_frag))))}))

  
  prepare=FALSE
}
#############################################################################################

pdf(file=paste(pathFigures, "eLife_Figures/Figure1_FigureSupplement1.pdf", sep=""), width=6.85, height=9)

m=matrix(rep(NA, 4*3), nrow=4)

m[1,]=1:3
m[2,]=4:6
m[3,]=7:9
m[4,]=10:12

layout(m)

#############################################################################################

labels=list("human"=c("a", "b","c"), "mouse"=c("d", "e", "f"))

## example plots for human and mouse

for(sp in c("human", "mouse")){
  sample=get(paste("exsample", sp, sep="."))
  obs=get(paste("obs", sp, sep="."))
  sim=get(paste("sim", sp, sep="."))

  this.info=get(paste("sampleinfo",sp,sep="."))
  
  celltype=this.info[which(this.info$Sample.ID==sample), "Broad.cell.type.or.tissue"]
  
  this.obs=obs[which(!is.na(obs[,sample])),]
  this.sim=sim[which(!is.na(sim[,sample])),]
  
  ## distance bait-contact
  
  d.obs=density(this.obs$distance, bw=0.5)
  d.sim=density(this.sim$distance, bw=0.5)
  
  xlim=range(c(0, maxDistance))
  ylim=range(c(d.obs$y, d.sim$y))
    
  xaxs=c(0, 500e3, 1e6, 1.5e6, 2e6, 2.5e6)
  xaxslabels=c("0", "0.5", "1", "1.5", "2", "2.5")  
  
  par(mar=c(3.1, 3.1, 2.1, 1.1))
  plot(d.obs$x, d.obs$y, type="l", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, col=dataset.colors["Original"])
  lines(d.sim$x, d.sim$y, col=dataset.colors["Simulated"])
  axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.95)
  mtext("distance bait-fragment (Mb)", side=1, line=1.5, cex=0.75)
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.95)
  
  mtext("density of nb contacts", side=2, line=2, cex=0.75)
  
  median.obs=round(median(this.obs$distance)/1000, digits=0)
  median.sim=round(median(this.sim$distance)/1000, digits=0)
  mtext(paste("median ", median.obs, " Kb", sep=""), side=3, line=-1.5, cex=0.75, col=dataset.colors["Original"])
  mtext(paste("median ", median.sim, " Kb", sep=""), side=3, line=-2.5, cex=0.75, col=dataset.colors["Simulated"])
  
  ## legend
  
  if(sp=="human"){
    legend("right", bty="n", col=dataset.colors, lty=1, legend=c("PCHi-C data", "simulated data"), cex=1.1, inset=0.05, xpd=NA)
  }

  ## label
  mtext(labels[[sp]][1], side=3, at=-0.5e6, cex=1.1, font=2, line=0.85)
  
  ## number of interactions per bait
  
  nb.contacts.per.bait.obs=as.numeric(table(this.obs$id_bait))
  nb.contacts.per.bait.sim=as.numeric(table(this.sim$id_bait))
  
  maxval=max(c(20, nb.contacts.per.bait.obs, nb.contacts.per.bait.sim))
  
  tab.obs=table(cut(nb.contacts.per.bait.obs, breaks=unique(sort(c(0:19, maxval))), include.lowest=T))
  tab.sim=table(cut(nb.contacts.per.bait.sim, breaks=unique(sort(c(0:19, maxval))), include.lowest=T))
  m=matrix(c(tab.obs, tab.sim), nrow=2, byrow=T)
  
  par(mar=c(3.1, 3.1, 2.1, 1.1))
  b=barplot(m, beside=T, col=dataset.colors, border=NA, axes=F, xlab="", ylab="")
  xpos=apply(b,2,mean)
  axis(side=1, at=xpos, labels=c(as.character(1:19), ""), cex=0.95, mgp=c(3, 0.5, 0))
  mtext("20+", at=xpos[length(xpos)]+diff(xpos)[1]/2, side=1, line=0.5, cex=0.65)
  
  axis(side=2, cex=0.95, mgp=c(3, 0.5, 0))
  mtext("nb. baits", side=2, line=1.5, cex=0.75)
  mtext("nb. contacted fragments", side=1, line=1.5, cex=0.75)
  
  mean.obs=round(mean(nb.contacts.per.bait.obs), digits=2)
  mtext(paste("mean ", mean.obs, sep=""), side=3, line=-1.5, cex=0.75, col=dataset.colors["Original"])
  mean.sim=round(mean(nb.contacts.per.bait.sim), digits=2)
  mtext(paste("mean ", mean.sim, sep=""), side=3, line=-2.5, cex=0.75, col=dataset.colors["Simulated"])

  ## labels
  mtext(labels[[sp]][2], side=3, at=-12, cex=1.1, font=2, line=0.85)
  
  ## number of interactions per fragment
  
  nb.contacts.per.frag.obs=as.numeric(table(this.obs$id_frag))
  nb.contacts.per.frag.sim=as.numeric(table(this.sim$id_frag))
  
  maxval=max(c(20, nb.contacts.per.frag.obs, nb.contacts.per.frag.sim))
  
  tab.obs=table(cut(nb.contacts.per.frag.obs, breaks=unique(sort(c(0:9, maxval))), include.lowest=T))
  tab.sim=table(cut(nb.contacts.per.frag.sim, breaks=unique(sort(c(0:9, maxval))), include.lowest=T))
  m=matrix(c(tab.obs, tab.sim), nrow=2, byrow=T)
  
  par(mar=c(3.1, 3.1, 2.1, 3.1))
  b=barplot(m, beside=T, col=dataset.colors, border=NA, axes=F, xlab="", ylab="")
  xpos=apply(b,2,mean)
  axis(side=1, at=xpos, labels=c(as.character(1:9), ""), cex=0.95, mgp=c(3, 0.5, 0))
  
  axis(side=2, cex=0.95, mgp=c(3, 0.5, 0))
  
  mtext("nb. fragments", side=2, line=1.5, cex=0.75)
  mtext("nb. contacting baits", side=1, line=1.5, cex=0.75)
  
  nb=dim(this.obs)[1]
  
  mtext(paste(sp,"\n",sample, "\n", "N=", nb,  sep=""), side=4, cex=0.75, las=2, line=-2) 
  
  mean.obs=round(mean(nb.contacts.per.frag.obs), digits=2)
  mtext(paste("mean ", mean.obs, sep=""), side=3, line=-1.5, cex=0.75, col=dataset.colors["Original"])
  mean.sim=round(mean(nb.contacts.per.frag.sim), digits=2)
  mtext(paste("mean ", mean.sim, sep=""), side=3, line=-2.5, cex=0.75, col=dataset.colors["Simulated"])

  mtext(labels[[sp]][3], side=3, at=-5.5, cex=1.1, font=2, line=0.85)
  
}

#############################################################################################

## median distance, human

par(mar=c(3.1, 3.1, 2.1, 1.1))

plot(median.distance.obs.human/1000, median.distance.sim.human/1000, pch=20, xlab="", ylab="", axes=F, main="", xlim=c(150, 400), ylim=c(150, 400))
abline(0,1)

axis(side=1, mgp=c(3, 0.5, 0))
axis(side=2, mgp=c(3, 0.5, 0))

box()

mtext("median dist. (Kb), PCHi-C data", side=1, line=1.75, cex=0.75)
mtext("median dist. (Kb), simulated", side=2, line=2.1, cex=0.75)

mtext("g", side=3, at=90, cex=1.1, font=2, line=1)

#############################################################################################

## mean nb contacts per bait

par(mar=c(3.1, 3.1, 2.1, 1.1))

plot(mean.nb.contacts.bait.obs.human, mean.nb.contacts.bait.sim.human, pch=20, xlab="", ylab="", axes=F, main="", xlim=c(3, 15), ylim=c(3, 15))
abline(0,1)

axis(side=1, mgp=c(3, 0.5, 0))
axis(side=2, mgp=c(3, 0.5, 0))

box()

mtext("mean nb. fragments, PCHi-C data", side=1, line=1.75, cex=0.75)
mtext("mean nb. fragments, simulated", side=2, line=2, cex=0.75)

mtext("h", side=3, at=0.1, cex=1.1, font=2, line=1)

#############################################################################################


## mean nb contacts per fragment

par(mar=c(3.1, 2.7, 2.1, 1.5))

plot(mean.nb.contacts.frag.obs.human, mean.nb.contacts.frag.sim.human, pch=20, xlab="", ylab="", axes=F, main="", xlim=c(1, 2), ylim=c(1,2))
abline(0,1)

axis(side=1, mgp=c(3, 0.5, 0))
axis(side=2, mgp=c(3, 0.5, 0))

box()

mtext("mean nb. baits, PCHi-C data", side=1, line=1.75, cex=0.75)
mtext("mean nb. baits, simulated", side=2, line=2, cex=0.75)

mtext("i", side=3, at=0.75, cex=1.1, font=2, line=1)

mtext("human", side=4, line=0.1, cex=0.75)

#############################################################################################
#############################################################################################

## median distance, mouse

par(mar=c(3.1, 3.1, 2.1, 1.1))

plot(median.distance.obs.mouse/1000, median.distance.sim.mouse/1000, pch=20, xlab="", ylab="", axes=F, main="")#, xlim=c(150, 400), ylim=c(150, 400))
abline(0,1)

axis(side=1, mgp=c(3, 0.5, 0))
axis(side=2, mgp=c(3, 0.5, 0))

box()

mtext("median dist. (Kb), PCHi-C data", side=1, line=1.75, cex=0.75)
mtext("median dist. (Kb), simulated", side=2, line=2.1, cex=0.75)

mtext("j", side=3, at=45, cex=1.1, font=2, line=1)

#############################################################################################

## mean nb contacts per bait

par(mar=c(3.1, 3.1, 2.1, 1.1))

plot(mean.nb.contacts.bait.obs.mouse, mean.nb.contacts.bait.sim.mouse, pch=20, xlab="", ylab="", axes=F, main="", xlim=c(4, 13), ylim=c(4, 13))
abline(0,1)

axis(side=1, mgp=c(3, 0.5, 0))
axis(side=2, mgp=c(3, 0.5, 0))

box()

mtext("mean nb. fragments, PCHi-C data", side=1, line=1.75, cex=0.75)
mtext("mean nb. fragments, simulated", side=2, line=2, cex=0.75)

mtext("k", side=3, at=1.9, cex=1.1, font=2, line=1)

#############################################################################################

## mean nb contacts per fragment

par(mar=c(3.1, 2.7, 2.1, 1.5))

plot(mean.nb.contacts.frag.obs.mouse, mean.nb.contacts.frag.sim.mouse, pch=20, xlab="", ylab="", axes=F, main="", xlim=c(1.1, 2.1), ylim=c(1.1,2.1))
abline(0,1)

axis(side=1, mgp=c(3, 0.5, 0))
axis(side=2, mgp=c(3, 0.5, 0))

box()

mtext("mean nb. baits, PCHi-C data", side=1, line=1.75, cex=0.75)
mtext("mean nb. baits, simulated", side=2, line=2, cex=0.75)

mtext("l", side=3, at=0.85, cex=1.1, font=2, line=1)

mtext("mouse", side=4, line=0.1, cex=0.75)

#############################################################################################


#############################################################################################

dev.off()

#############################################################################################
