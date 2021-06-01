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
  load(paste(pathFigures, "RData/data.sample.info.RData",sep=""))

   load=F
 }

#################################################################################

if(prepare){

  nb.unique.fragments.obs=list()
  nb.unique.fragments.sim=list()

  length.unique.fragments.obs=list()
  length.unique.fragments.sim=list()
  
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
    
  }
  
  prepare=F
}

#################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#############################################################################

pdf(paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_S3.pdf", sep=""), width=6.85, height=6.5)

m=matrix(1:4, nrow=2, byrow=T)
layout(m)

#############################################################################

## nb of contacted fragments, obs vs. sim

labels=c("A","B")
names(labels)=c("human", "mouse")

for(sp in c("human", "mouse")){

  par(mar=c(3.5, 4.1, 2.25, 1.5))

  lim=range(c(nb.unique.fragments.obs[[sp]], nb.unique.fragments.sim[[sp]]))/1000
  lim=lim+c(-diff(lim)/20, diff(lim)/20)
  
  plot(nb.unique.fragments.obs[[sp]]/1000, nb.unique.fragments.sim[[sp]]/1000, pch=20, xlab="", ylab="", axes=F, xlim=lim, ylim=lim)

  box()

  abline(0,1, lty=3)

  ax=pretty(lim)
  axis(side=1, mgp=c(3, 0.5, 0), cex.axis=0.95, at=ax, labels=paste(ax,"k",sep=""))
  mtext("nb. contacted fragments, PCHi-C", side=1, cex=0.85, line=1.75)

  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.95, at=ax, labels=paste(ax,"k",sep=""))
  mtext("nb. contacted fragments, simulations", side=2, cex=0.85, line=2.25)

  mtext(labels[sp], side=3, at=lim[1]-diff(lim)/4.48, font=2, line=1.15, cex=1.1)

  mtext(sp, side=3, line=0.5, cex=0.9)
}

#############################################################################

## length of contacted fragments, obs vs. sim

labels=c("C","D")
names(labels)=c("human", "mouse")

for(sp in c("human", "mouse")){

  par(mar=c(3.5, 4.1, 2.25, 1.5))

  lim=range(c(length.unique.fragments.obs[[sp]], length.unique.fragments.sim[[sp]]))/1e6

  lim=lim+c(-diff(lim)/20, diff(lim)/20)
  
  plot(length.unique.fragments.obs[[sp]]/1e6, length.unique.fragments.sim[[sp]]/1e6, pch=20, xlab="", ylab="", axes=F, xlim=lim, ylim=lim)

  box()

  abline(0,1, lty=3)

  ax=pretty(lim)
  axis(side=1, mgp=c(3, 0.5, 0), cex.axis=0.95, at=ax, labels=paste(ax,"Mb",sep=""))
  mtext("total length fragments, PCHi-C", side=1, cex=0.85, line=1.75)

  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.95, at=ax, labels=paste(ax,"Mb",sep=""))
  mtext("total length fragments, simulations", side=2, cex=0.85, line=2.25)

  mtext(labels[sp], side=3, at=lim[1]-diff(lim)/4.48, font=2, line=1.15, cex=1.1)

  mtext(sp, side=3, line=0.5, cex=0.9)

}

#############################################################################

dev.off()

#############################################################################
