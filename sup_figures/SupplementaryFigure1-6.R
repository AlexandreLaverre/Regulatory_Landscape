#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

##########################################################################

if(load){
 
  load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

  load=FALSE
}

##########################################################################

if(prepare){

  ## determine which samples we show in which figure
  ## we don't mix the two species, 7 samples per figure
  
  nbsamples.human=dim(sampleinfo[["human"]])[1]
  nbsamples.mouse=dim(sampleinfo[["mouse"]])[1]
  nbsamples=nbsamples.human+nbsamples.mouse

  nbperfile=7
  nbfiles.human=ceiling(nbsamples.human/nbperfile)
  nbfiles.mouse=ceiling(nbsamples.mouse/nbperfile)
  
  all.info=rbind(sampleinfo[["human"]], sampleinfo[["mouse"]], stringsAsFactors=F)

  all.info$PlotIndex=rep(NA, dim(all.info)[1])
  all.info$PlotIndex[which(all.info$Species=="human")]=rep(1:nbfiles.human, each=nbperfile)[1:nbsamples.human]
  all.info$PlotIndex[which(all.info$Species=="mouse")]=rep(nbfiles.human+(1:nbfiles.mouse), each=nbperfile)[1:nbsamples.mouse]

  nbfigures=max(all.info$PlotIndex)
  
  prepare=F
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#############################################################################

for(i in 1:nbfigures){
  print(i)

  this.info=all.info[which(all.info$PlotIndex==i),]
  this.sp=unique(this.info$Species)
  
  if(length(this.sp)!=1){
    stop("Weird! multiple species in this figure")
  }

  labels=toupper(letters[1:nbperfile])
  names(labels)=this.info$Sample.ID
  
  ## get observed and simulated contacts for this species
  obs=observed.contacts[[this.sp]]
  sim=simulated.contacts[[this.sp]]

  ## start figure
  pdf(paste(pathFigures, "SupplementaryFigure",i,".pdf", sep=""), width=6.85, height=11)
  
  ## layout
  
  ncol=3
  nrow=nbperfile
  m=matrix(1:(nrow*ncol), nrow=nrow, byrow=T)
  layout(m)
  
  par(oma=c(1, 0, 0, 0))
    
  for(sample in this.info$Sample.ID){
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
    
    par(mar=c(2.1, 3.1, 2.1, 1.1))
    plot(d.obs$x, d.obs$y, type="l", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, col=dataset.colors["Original"])
    lines(d.sim$x, d.sim$y, col=dataset.colors["Simulated"])
    axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.95)
    mtext("distance bait-fragment (Mb)", side=1, line=1.5, cex=0.7)
    
    axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.95)

    mtext("density of nb contacts", side=2, line=2, cex=0.7)

    median.obs=round(median(this.obs$distance)/1000, digits=0)
    median.sim=round(median(this.sim$distance)/1000, digits=0)
    mtext(paste("median ", median.obs, " Kb", sep=""), side=3, line=-1.5, cex=0.8, col=dataset.colors["Original"])
    mtext(paste("median ", median.sim, " Kb", sep=""), side=3, line=-2.5, cex=0.8, col=dataset.colors["Simulated"])
    
    ## legend

    if(sample==this.info$Sample.ID[1]){
      legend("right", bty="n", col=dataset.colors, lty=1, legend=c("PC-HiC data", "simulated data"), cex=1.1, inset=0.05, xpd=NA)
    }

    mtext(labels[sample], side=3, at=-0.55e6, cex=1.1, font=2, line=0.85)

    ## number of interactions per bait

    nb.contacts.per.bait.obs=as.numeric(table(this.obs$id_bait))
    nb.contacts.per.bait.sim=as.numeric(table(this.sim$id_bait))

    maxval=max(c(20, nb.contacts.per.bait.obs, nb.contacts.per.bait.sim))

    tab.obs=table(cut(nb.contacts.per.bait.obs, breaks=unique(sort(c(0:19, maxval))), include.lowest=T))
    tab.sim=table(cut(nb.contacts.per.bait.sim, breaks=unique(sort(c(0:19, maxval))), include.lowest=T))
    m=matrix(c(tab.obs, tab.sim), nrow=2, byrow=T)

    par(mar=c(2.1, 3.1, 2.1, 1.1))
    b=barplot(m, beside=T, col=dataset.colors, border=NA, axes=F, xlab="", ylab="")
    xpos=apply(b,2,mean)
    axis(side=1, at=xpos, labels=c(as.character(1:19), ""), cex=0.95, mgp=c(3, 0.5, 0))
    mtext("20+", at=xpos[length(xpos)]+diff(xpos)[1]/2, side=1, line=0.5, cex=0.6)
    
    axis(side=2, cex=0.95, mgp=c(3, 0.5, 0))
    mtext("nb. baits", side=2, line=1.5, cex=0.7)
    mtext("nb. contacted fragments", side=1, line=1.5, cex=0.7)

    median.obs=median(nb.contacts.per.bait.obs)
    mtext(paste("median ", median.obs, sep=""), side=3, line=-1.5, cex=0.8, col=dataset.colors["Original"])

  
    ## number of interactions per fragment

    nb.contacts.per.frag.obs=as.numeric(table(this.obs$id_frag))
    nb.contacts.per.frag.sim=as.numeric(table(this.sim$id_frag))

    maxval=max(c(20, nb.contacts.per.frag.obs, nb.contacts.per.frag.sim))

    tab.obs=table(cut(nb.contacts.per.frag.obs, breaks=unique(sort(c(0:9, maxval))), include.lowest=T))
    tab.sim=table(cut(nb.contacts.per.frag.sim, breaks=unique(sort(c(0:9, maxval))), include.lowest=T))
    m=matrix(c(tab.obs, tab.sim), nrow=2, byrow=T)

    par(mar=c(2.1, 3.1, 2.1, 4.1))
    b=barplot(m, beside=T, col=dataset.colors, border=NA, axes=F, xlab="", ylab="")
    xpos=apply(b,2,mean)
    axis(side=1, at=xpos, labels=c(as.character(1:9), ""), cex=0.95, mgp=c(3, 0.5, 0))
   
    axis(side=2, cex=0.95, mgp=c(3, 0.5, 0))
   
    mtext("nb. fragments", side=2, line=1.5, cex=0.7)
    mtext("nb. contacting baits", side=1, line=1.5, cex=0.7)

    nb=dim(this.obs)[1]

    mtext(paste(this.sp,"\n",sample, "\n", "N=", nb,  sep=""), side=4, cex=0.75, las=2, line=-2) 

    median.obs=median(nb.contacts.per.frag.obs)
    mtext(paste("median ", median.obs, sep=""), side=3, line=-1.5, cex=0.8, col=dataset.colors["Original"])
    median.sim=median(nb.contacts.per.frag.sim)
    mtext(paste("median ", median.obs, sep=""), side=3, line=-2.5, cex=0.8, col=dataset.colors["Simulated"])

  }

  dev.off()

}

##########################################################################

