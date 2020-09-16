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
  this.info=all.info[which(all.info$PlotIndex==i),]
  this.sp=unique(this.info$Species)

  if(length(this.sp)!=1){
    stop("Weird! multiple species in this figure")
  }

  ## get observed and simulated contacts for this species
  obs=observed.contacts[[this.sp]]
  sim=simulated.contacts[[this.sp]]


  ## start figure
  pdf(paste(pathFigures, "SupplementaryFigure",i,".pdf", sep=""), width=6.85, height=11)
  
  ## layout
  
  ncol=3
  nrow=nbperfile
  m=matrix(1:(nrow*ncol), nrow=nrow)
  layout(m)
  
  par(oma=c(2, 2, 1, 1))
    
  for(sample in this.info$Sample.ID){
    this.obs=obs[which(!is.na(obs[,sample])),]
    this.sim=sim[which(!is.na(sim[,sample])),]

    ## distance bait-contact
    
    d.obs=density(this.obs$distance, bw=0.1)
    d.sim=density(this.sim$distance, bw=0.1)
    
    xlim=range(c(0, maxDistance))
    ylim=range(c(d.obs$y, d.sim$y))
    
    xaxs=c(0, 500e3, 1e6, 1.5e6, 2e6, 2.5e6)
    xaxslabels=c("0", "0.5", "1", "1.5M", "2", "2.5")  
    
    par(mar=c(2.1, 2.1, 1.1, 1.1))
    plot(d.obs$x, d.obs$y, type="l", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, col=dataset.colors["Original"])
    lines(d.sim$x, d.sim$y, col=dataset.colors["Simulated"])
    axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.6)
    axis(side=2, mgp=c(3, 0.5, 0), cex.axis=0.6)

    ## number of interactions per bait

    nb.contacts.per.bait.obs=as.numeric(table(this.obs$id_bait))
    nb.contacts.per.bait.sim=as.numeric(table(this.sim$id_bait))

    d.obs=density(nb.contacts.per.bait.obs, bw=0.1)
    d.sim=density(nb.contacts.per.bait.sim, bw=0.1)
    
    xlim=range(c(0, max(c(nb.contacts.per.bait.obs, nb.contacts.per.bait.sim))))
    ylim=range(c(d.obs$y, d.sim$y))
    
    xaxs=c(0, 500e3, 1e6, 1.5e6, 2e6, 2.5e6)
    xaxslabels=c("0", "0.5", "1", "1.5M", "2", "2.5")  
    
    par(mar=c(2.1, 2.1, 1.1, 1.1))
    plot(d.obs$x, d.obs$y, type="l", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, col=dataset.colors["Original"])
    lines(d.sim$x, d.sim$y, col=dataset.colors["Simulated"])
    axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.6)
    axis(side=2, mgp=c(3, 0.5, 0), cex.axis=0.6)
    
    ## number of interactions per fragment

    nb.contacts.per.frag.obs=as.numeric(table(this.obs$id_frag))
    nb.contacts.per.frag.sim=as.numeric(table(this.sim$id_frag))

    d.obs=density(nb.contacts.per.frag.obs, bw=0.1)
    d.sim=density(nb.contacts.per.frag.sim, bw=0.1)
    
    xlim=range(c(0, max(c(nb.contacts.per.bait.obs, nb.contacts.per.bait.sim))))
    ylim=range(c(d.obs$y, d.sim$y))
    
    xaxs=c(0, 500e3, 1e6, 1.5e6, 2e6, 2.5e6)
    xaxslabels=c("0", "0.5", "1", "1.5M", "2", "2.5")  
    
    par(mar=c(2.1, 2.1, 1.1, 1.1))
    plot(d.obs$x, d.obs$y, type="l", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, col=dataset.colors["Original"])
    lines(d.sim$x, d.sim$y, col=dataset.colors["Simulated"])
    axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.6)
    axis(side=2, mgp=c(3, 0.5, 0), cex.axis=0.6)
    
    
  }


  dev.off()

}

##########################################################################

