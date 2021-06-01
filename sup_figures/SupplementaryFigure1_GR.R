#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=TRUE
  prepare=TRUE
  source("../main_figures/parameters.R")

  pathEnhancers="../../../RegulatoryLandscapesManuscript/SupplementaryDataset4/"

  library(plotrix)
}

##########################################################################

if(load){
  sp="human"

  load(paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))

  obs=observed.contacts[[sp]]
  sim=simulated.contacts[[sp]]

  load(paste(pathFigures, "RData/data.bait.annotation.RData",sep=""))

  baits=bait.info[[sp]]

  load(paste(pathFigures, "RData/data.restriction.map.RData", sep=""))

  map=restriction.map[[sp]]

  ## enhancer coordinates, ENCODE
  enhancers=read.table(paste(pathEnhancers, sp, "/ENCODE/enhancer_coordinates.bed", sep=""), h=T, stringsAsFactors=F)
  
  load=FALSE
}

##########################################################################

if(prepare){
  ## "chr1:147541410:147546578"
  ## "chrX:30722815:30725270"
  
  ## bait="chr4:41081621:41085450"
  ## cell="PEK_early"
  
  ## bait="chr12:53218969:53240775"
  ## cell="hESC"

  bait="chr11:66263828:66266552"
  cell="CD34"

  this.obs=obs[which(obs$id_bait==bait & !is.na(obs[,cell])),]
  this.sim=sim[which(sim$id_bait==bait & !is.na(sim[,cell])),]

  chr=this.obs$chr_bait[1]

  xrange=range(c(this.obs$start, this.sim$start, this.obs$end, this.sim$end))

  this.map=map[which(map$chr==chr & map$start>=xrange[1] & map$end<=xrange[2]),]

  this.baits=baits[which(baits$chr==chr & baits$start>=xrange[1] & baits$end<=xrange[2]),]

  this.enhancers=enhancers[which(enhancers$chr==chr & enhancers$start>=xrange[1] & enhancers$end<=xrange[2]),]


  ## all obs, only this sample

  all.obs=obs[which(!is.na(obs[,cell])),]
  all.sim=sim[which(!is.na(sim[,cell])),]
  
  prepare=FALSE
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_S1.pdf", sep=""), width=6.85, height=3)

m=matrix(rep(NA, 10*40), nrow=10)

for(i in 1){
  m[i,]=c(rep(1, 25), rep(6, 15))
}

for(i in 2:4){
  m[i,]=c(rep(2, 25), rep(6, 15))
}

for(i in 5){
 m[i,]=c(rep(3, 25), rep(6, 15))
}

for(i in 6:9){
 m[i,]=c(rep(4, 25), rep(7, 15))
}

for(i in 10){
 m[i,]=c(rep(5, 25), rep(7, 15))
}

layout(m)

##########################################################################
## margin

par(mar=c(0.2, 5.2, 0.2, 1.2))

##########################################################################

## legend plot

plot(1, type="n", xlim=c(0,1.15), ylim=c(0,1), axes=F, xlab="", ylab="", xaxs="i", yaxs="i")

legend("topright", c("PCHi-C data", "simulated data"), col=dataset.colors, lty=1, inset=0.01, bty="n", xpd=NA, cex=0.95, seg.len=1)

small=0.01

rect(0.815, -0.75, 0.83, -0.55, col="gray20", border="gray20", xpd=NA)
rect(0.84, -0.75, 0.85, -0.55, col="gray80", border="gray80", xpd=NA)

text(x=0.875, y=-0.6, labels="restriction fragments", cex=0.95, xpd=NA, adj=c(0,0.5))

rect(0.83, -1.25, 0.832, -1.05, col="gray20", border="gray20", xpd=NA)
rect(0.84, -1.25, 0.842, -1.05, col="gray20", border="gray20", xpd=NA)
rect(0.82, -1.25, 0.822, -1.05, col="gray20", border="gray20", xpd=NA)

text(x=0.875, y=-1.125, labels="enhancers", cex=0.95, xpd=NA, adj=c(0,0.5))

rect(0.83, -1.75, 0.845, -1.55, col="gray80", border="red", xpd=NA)

text(x=0.875, y=-1.65, labels="bait", cex=0.95, xpd=NA, adj=c(0,0.5))

mtext("example bait, one sample", side=3, line=-1, at=0.05, cex=0.65, font=2)

##########################################################################

## observed interactions

yrange=c(0, 1)

plot(1, type="n", xlim=xrange, ylim=c(0,0.6), axes=F, xlab="", ylab="", xaxs="i", yaxs="i")

for(i in 1:nrow(this.obs)){
  pos.bait=(this.obs$start_bait[i]+this.obs$end_bait[i])/2
  pos.frag=(this.obs$start[i]+this.obs$end[i])/2
  x.center=(pos.bait+pos.frag)/2
  R.x=abs(x.center-pos.bait)
  R.y=0.5
  
  draw.ellipse(x=x.center, y=0, a=R.x, b=0.5, col=NA, border=dataset.colors["Original"])
}

##########################################################################

## restriction fragments

plot(1, type="n", xlim=xrange, ylim=c(-1,1), axes=F, xlab="", ylab="", xaxs="i")

for(i in 1:dim(this.map)[1]){
  if(i%%2==0){
    rect(this.map$start[i], 0.25, this.map$end[i], 0.75, col="gray80", border=NA)
  } else{
    rect(this.map$start[i], 0.25, this.map$end[i], 0.75, col="gray20", border=NA)
  }
}

mtext("fragments", side=2, las=2, at=0.55, cex=0.65,line=1)

## coordinates for the bait

rect(this.obs$start_bait[1], 0.25, this.obs$end_bait[1], 0.75, col="gray80", border="red")
abline(h=-0.5, lty=1, col="gray80")

for(i in 1:dim(this.enhancers)[1]){
  rect(this.enhancers$start[i], -0.25, this.enhancers$end[i], -0.75, col="gray20", border=NA)
}

mtext("enhancers", side=2, las=2, at=-0.45, cex=0.65,line=1)

##########################################################################

## simulated interactions

yrange=c(0, 1)

plot(1, type="n", xlim=xrange, ylim=c(-0.6,0), axes=F, xlab="", ylab="", xaxs="i", yaxs="i")

for(i in 1:nrow(this.sim)){
  pos.bait=(this.sim$start_bait[i]+this.sim$end_bait[i])/2
  pos.frag=(this.sim$start[i]+this.sim$end[i])/2
  x.center=(pos.bait+pos.frag)/2
  R.x=abs(x.center-pos.bait)
  R.y=0.5
  
  draw.ellipse(x=x.center, y=0, a=R.x, b=-0.5, col=NA, border=dataset.colors["Simulated"])
  
}

##########################################################################

plot(1, type="n", xlim=xrange, ylim=c(-0.6,0), axes=F, xlab="", ylab="", xaxs="i", yaxs="i")

xaxs=pretty(xrange/1e6, n=10)
axis(side=1, mgp=c(3, 0.5, 0), cex.axis=0.9, line=-1.5, at=xaxs*1e6, labels=paste0(xaxs, "Mb"))

mtext(chr, side=1, line=-2, at=xrange[1]-diff(xrange)/20, cex=0.625)

##########################################################################

## density plot, distance, observed interactions

d.obs=density(all.obs$distance, bw=0.5)
d.sim=density(all.sim$distance, bw=0.5)

xlim=range(c(0, maxDistance))
ylim=range(c(d.obs$y, d.sim$y))

xaxs=c(0, 500e3, 1e6, 1.5e6, 2e6, 2.5e6)
xaxslabels=c("0", "0.5", "1", "1.5", "2", "2.5")  

par(mar=c(2.8, 7.1, 1.1, 0.25))
plot(d.obs$x, d.obs$y, type="l", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, col=dataset.colors["Original"])
axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.9)
mtext("distance bait-fragment (Mb)", side=1, line=1.5, cex=0.65)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.9)

mtext("density, nb. contacts", side=2, line=2.1, cex=0.65)

mtext("all baits, one sample", side=3, line=0, at=xlim[2], adj=1, cex=0.65, font=2)

##########################################################################

## density plot, distance, simulated interactions

plot(d.sim$x, d.sim$y, type="l", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, col=dataset.colors["Simulated"])
axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.9)
mtext("distance bait-fragment (Mb)", side=1, line=1.5, cex=0.65)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.9)

mtext("density, nb. contacts", side=2, line=2.1, cex=0.65)

dev.off()

##########################################################################
