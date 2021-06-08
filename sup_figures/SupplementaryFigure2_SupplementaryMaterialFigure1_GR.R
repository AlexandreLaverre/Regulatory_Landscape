#################################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

#################################################################################

if(load){

  load(paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))
  load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))
  load(paste(pathFigures, "RData/data.sample.info.RData",sep=""))

   load=F
 }

#################################################################################

if(prepare){

  nb.unique.fragments.obs=list()
  nb.unique.fragments.sim=list()

  length.unique.fragments.obs=list()
  length.unique.fragments.sim=list()

  pc.degree.frag=list()
  pc.degree.bait=list()
  
  for(sp in c("human", "mouse")){
    print(sp)
    
    samples=sampleinfo[[sp]][,"Sample.ID"]
  
    obs=observed.contacts[[sp]]
    sim=simulated.contacts[[sp]]
    
    nb.unique.fragments.obs[[sp]]=unlist(lapply(samples, function(x) length(unique(obs$id_frag[which(obs[,x]>0)]))))
    nb.unique.fragments.sim[[sp]]=unlist(lapply(samples, function(x) length(unique(sim$id_frag[which(obs[,x]>0)]))))
    
    length.unique.fragments.obs[[sp]]=unlist(lapply(samples, function(x) {y=obs[which(obs[,x]>0), ]; y=y[which(!duplicated(y$id_frag)),]; return(sum(y$contacted_length))}))
    length.unique.fragments.sim[[sp]]=unlist(lapply(samples, function(x) {y=sim[which(sim[,x]>0), ]; y=y[which(!duplicated(y$id_frag)),]; return(sum(y$contacted_length))}))

    print(paste(sum(obs$contacted_length[which(!duplicated(obs$id_frag))])/1e6, "Mb total covered length, observed data"))
    print(paste(sum(sim$contacted_length[which(!duplicated(sim$id_frag))])/1e6, "Mb total covered length, simulated data"))

    nb.obs=length(unique(obs$id_frag))
    nb.sim=length(unique(sim$id_frag))
    nb.common=length(intersect(obs$id_frag, sim$id_frag))

    print(paste("observed",nb.obs,"fragments, simulated",nb.sim,"fragments,", nb.common, "in common"))

    obs$id=paste(obs$id_bait, obs$id_frag, sep="-")
    sim$id=paste(sim$id_bait, sim$id_frag, sep="-")

    print(paste("observed",nrow(obs),"contacts, simulated",nrow(sim),"contacts,", length(intersect(obs$id, sim$id)), "in common"))


    
    stats.obs=fragment.statistics[[sp]][["original"]]
    stats.sim=fragment.statistics[[sp]][["simulated"]]
    
  ## number of contacts per bait in real and simulated data, after filtering
  
    bait.degree.obs=table(as.factor(obs$id_bait))
    bait.degree.sim=table(as.factor(sim$id_bait))
    
    tab.bait.degree.obs=table(cut(bait.degree.obs, breaks=c(seq(from=0, to=50, by=5), max(bait.degree.obs)), include.lowest=T))
    tab.bait.degree.sim=table(cut(bait.degree.sim, breaks=c(seq(from=0, to=50, by=5), max(bait.degree.sim)), include.lowest=T))
    
    this.pc.degree.bait=matrix(c(tab.bait.degree.obs, tab.bait.degree.sim), nrow=2, byrow=T)
    this.pc.degree.bait=100*this.pc.degree.bait/apply(this.pc.degree.bait,1,sum)
    
    colnames(this.pc.degree.bait)=paste(seq(from=1, to=50, by=5), seq(from=5, to=55, by=5), sep="-")
    colnames(this.pc.degree.bait)[ncol(this.pc.degree.bait)]=">50"
    
    ## number of contacts per fragment in real and simulated data, after filtering
    
    frag.degree.obs=table(as.factor(obs$id_frag))
    frag.degree.sim=table(as.factor(sim$id_frag))
    
    tab.frag.degree.obs=table(cut(frag.degree.obs, breaks=c(seq(from=0, to=10, by=1), max(frag.degree.obs)), include.lowest=T))
    tab.frag.degree.sim=table(cut(frag.degree.sim, breaks=c(seq(from=0, to=10, by=1), max(frag.degree.sim)), include.lowest=T))
    
    this.pc.degree.frag=matrix(c(tab.frag.degree.obs, tab.frag.degree.sim), nrow=2, byrow=T)
    this.pc.degree.frag=100*this.pc.degree.frag/apply(this.pc.degree.frag,1,sum)
    
    colnames(this.pc.degree.frag)=as.character(1:11)
    colnames(this.pc.degree.frag)[ncol(this.pc.degree.frag)]=">10"

    pc.degree.bait[[sp]]=this.pc.degree.bait
    pc.degree.frag[[sp]]=this.pc.degree.frag
    
  }
  
  prepare=F
}

#################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#############################################################################

for(sp in c("human", "mouse")){

  if(sp=="human"){
    pdf(paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_S3.pdf", sep=""), width=5.5, height=5.75)
  } else{
    pdf(paste(pathFigures, "GenomeResearch_Figures/SupplementaryMaterialFigure1.pdf", sep=""), width=5.5, height=5.75)
  }
   
  m=matrix(1:4, nrow=2, byrow=T)
  layout(m)
  
#############################################################################
  
  ## nb of contacted fragments, obs vs. sim
  
 
  par(mar=c(4.75, 2.85, 2.1, 1.5))

  lim=range(c(nb.unique.fragments.obs[[sp]], nb.unique.fragments.sim[[sp]]))/1000
  lim=lim+c(-diff(lim)/20, diff(lim)/20)
  
  plot(nb.unique.fragments.obs[[sp]]/1000, nb.unique.fragments.sim[[sp]]/1000, pch=20, xlab="", ylab="", axes=F, xlim=lim, ylim=lim)

  box()

  abline(0,1, lty=3)

  ax=pretty(lim)
  axis(side=1, mgp=c(3, 0.5, 0), cex.axis=0.8, at=ax, labels=ax)
  mtext("nb. K contacted frag., PCHi-C", side=1, cex=0.75, line=1.75)

  axis(side=2, mgp=c(3, 0.5, 0), cex.axis=0.8, at=ax, labels=ax)
  mtext("nb. K contacted frag., simulations", side=2, cex=0.75, line=1.75)

  mtext("A", side=3, at=lim[1]-diff(lim)/4.35, font=2, line=1, cex=1)

#############################################################################

## length of contacted fragments, obs vs. sim
  
  par(mar=c(4.75, 2.85, 2.1, 1.5))
  
  lim=range(c(length.unique.fragments.obs[[sp]], length.unique.fragments.sim[[sp]]))/1e6

  lim=lim+c(-diff(lim)/20, diff(lim)/20)
  
  plot(length.unique.fragments.obs[[sp]]/1e6, length.unique.fragments.sim[[sp]]/1e6, pch=20, xlab="", ylab="", axes=F, xlim=lim, ylim=lim)

  box()

  abline(0,1, lty=3)

  ax=pretty(lim)
  axis(side=1, mgp=c(3, 0.5, 0), cex.axis=0.8, at=ax, labels=ax)
  mtext("tot. length frag., PCHi-C (Mb)", side=1, cex=0.75, line=1.75)

  axis(side=2, mgp=c(3, 0.5, 0), cex.axis=0.8, at=ax, labels=ax)
  mtext("tot. length frag., simulations (Mb)", side=2, cex=0.75, line=1.75)

  mtext("B", side=3, at=lim[1]-diff(lim)/4.35, font=2, line=1, cex=1.1)
  
#############################################################################
  
  ## bait degree, all samples
  
  par(mar=c(3.5, 2.75, 1.5, 0.5))
  
  ylim=c(0, max(as.numeric(pc.degree.bait[[sp]])))
  
  b=barplot(pc.degree.bait[[sp]], beside=T, xlab='', ylim=ylim, space=c(0.4,1), names=rep("", dim(pc.degree.bait[[sp]])[2]), ylab="", border=dataset.colors[c("Original", "Simulated")], col=dataset.colors[c("Original", "Simulated")], lwd=1.5,  mgp=c(3, 0.75, 0), axes=F, las=1)

  axis(side=2, mgp=c(3, 0.5, 0), cex.axis=0.8)
  
  ## axis labels
  mtext("number of contacted fragments", side=1, line=2.25, cex=0.75)
  mtext("% baits", side=2, line=1.5, cex=0.75)
  
  mtext(colnames(pc.degree.bait[[sp]]), at=apply(b, 2, mean), line=0.25, cex=0.7, side=1, las=2)
  
  ## legend & plot label
  legend("topleft", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],  fill=dataset.colors[c("Original", "Simulated")], bty='n', cex=0.85, inset=c(0.05, -0.15), xpd=NA)
  
  mtext("C", side=3, line=1, at=-6, font=2, cex=1.1)

##########################################################################
  
  ## fragment degree, all samples
  
  par(mar=c(3.5, 2.75, 1.5, 0.5))

  ylim=c(0, max(as.numeric(pc.degree.frag[[sp]])))
  
  b=barplot(pc.degree.frag[[sp]], beside=T, xlab='', ylim=ylim, space=c(0.4,1), names=rep("", dim(pc.degree.frag[[sp]])[2]), ylab="", border=dataset.colors[c("Original", "Simulated")], col=dataset.colors[c("Original", "Simulated")], lwd=1.5,  mgp=c(3, 0.75, 0), axes=F)

   axis(side=2, mgp=c(3, 0.5, 0), cex.axis=0.8)
  
  ## axis labels
  mtext("number of contacting baits", side=1, line=1.75, cex=0.75)
  mtext("% fragments", side=2, line=1.5, cex=0.75)
  
  mtext(colnames(pc.degree.frag[[sp]]), at=apply(b, 2, mean), line=0.25, cex=0.7, side=1, las=1)
  
  ## plot label
  
  mtext("D", side=3, line=1, at=-6, font=2, cex=1.1)
  
#############################################################################
  
  dev.off()

#############################################################################

}
