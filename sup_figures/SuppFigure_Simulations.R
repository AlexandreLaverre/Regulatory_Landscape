#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=TRUE
  prepare=TRUE
  source("../main_figures/parameters.R")

  library(plotrix)
  library(DescTools)
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
    
  prepare=FALSE
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "SupplementaryFigure1_Simulations.pdf", sep=""), width=6.85, height=2.5)

m=matrix(rep(NA, 9*40), nrow=9)

for(i in 1:3){
  m[i,]=rep(1, 40)
}

for(i in 4){
 m[i,]=rep(2, 40)
}

for(i in 5:8){
 m[i,]=rep(3, 40)
}

for(i in 9){
 m[i,]=rep(4, 40)
}

layout(m)

##########################################################################

## observed interactions

par(mar=c(0.2, 1.2, 0.2, 1.2))

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

plot(1, type="n", xlim=xrange, ylim=c(0,1), axes=F, xlab="", ylab="", xaxs="i")

for(i in 1:dim(this.map)[1]){
  if(i%%2==0){
    rect(this.map$start[i], 0.25, this.map$end[i], 0.75, col="gray80", border=NA)
  } else{
    rect(this.map$start[i], 0.25, this.map$end[i], 0.75, col="gray20", border=NA)
  }
}

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
par(mar=c(0.4, 1.2, 0.0, 1.2))

xaxs=pretty(xrange/1e6, n=10)
#xaxs=xaxs[which(xaxs>=(xrange[1]/1e6) & xaxs<=(xrange[2]/1e6))]
axis(side=1, mgp=c(3, 0.5, 0), cex=0.5, line=-1.5, at=xaxs*1e6, labels=paste0(xaxs, "Mb"))

mtext(chr, side=1, line=-2.75, at=xrange[1]+diff(xrange)/40, cex=0.75)

##########################################################################

dev.off()

##########################################################################
