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
  sp="human"
  
  load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))

  contact.obs=observed.contacts[[sp]]
  contact.sim=simulated.contacts[[sp]]

  stats.obs=fragment.statistics[[sp]][["original"]]
  stats.sim=fragment.statistics[[sp]][["simulated"]]

  load=FALSE
}

##########################################################################

if(prepare){
  ## number of contacts per bait in real and simulated data, after filtering
  
  bait.degree.obs=table(as.factor(contact.obs$id_bait))
  bait.degree.sim=table(as.factor(contact.sim$id_bait))

  tab.bait.degree.obs=table(cut(bait.degree.obs, breaks=c(seq(from=0, to=50, by=5), max(bait.degree.obs)), include.lowest=T))
  tab.bait.degree.sim=table(cut(bait.degree.sim, breaks=c(seq(from=0, to=50, by=5), max(bait.degree.sim)), include.lowest=T))

  pc.degree.bait=matrix(c(tab.bait.degree.obs, tab.bait.degree.sim), nrow=2, byrow=T)
  pc.degree.bait=100*pc.degree.bait/apply(pc.degree.bait,1,sum)

  colnames(pc.degree.bait)=paste(seq(from=1, to=50, by=5), seq(from=5, to=55, by=5), sep="-")
  colnames(pc.degree.bait)[ncol(pc.degree.bait)]=">50"

  ## number of contacts per fragment in real and simulated data, after filtering
  
  frag.degree.obs=table(as.factor(contact.obs$id_frag))
  frag.degree.sim=table(as.factor(contact.sim$id_frag))
 
  tab.frag.degree.obs=table(cut(frag.degree.obs, breaks=c(seq(from=0, to=10, by=1), max(frag.degree.obs)), include.lowest=T))
  tab.frag.degree.sim=table(cut(frag.degree.sim, breaks=c(seq(from=0, to=10, by=1), max(frag.degree.sim)), include.lowest=T))

  pc.degree.frag=matrix(c(tab.frag.degree.obs, tab.frag.degree.sim), nrow=2, byrow=T)
  pc.degree.frag=100*pc.degree.frag/apply(pc.degree.frag,1,sum)

  colnames(pc.degree.frag)=as.character(1:11)
  colnames(pc.degree.frag)[ncol(pc.degree.frag)]=">10"

   
  prepare=FALSE
}

##########################################################################

## this figure shows some statistics for observed and simulated interactions
## human only
## bait degree
## fragment degree
## correlation between nb of observed contacts and nb of genes/baits in the neighborhood of the fragments, observed & simulated
## gene density for contacted fragments and simulated fragments, observed & simulated

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "SupplementaryMaterialFigure8.pdf", sep=""), width=6.85, height=6.5)

m=matrix(1:4, nrow=2, byrow=T)

layout(m)

##########################################################################

## bait degree

par(mar=c(4.5, 4.5, 2.1, 1.1))

ylim=c(0, max(as.numeric(pc.degree.bait)))

b=barplot(pc.degree.bait, beside=T, xlab='', ylim=ylim, space=c(0.4,1), names=rep("", dim(pc.degree.bait)[2]), ylab="", border=dataset.colors[c("Original", "Simulated")], col=dataset.colors[c("Original", "Simulated")], lwd=1.5,  mgp=c(3, 0.75, 0), cex.axis=1, las=2)

## axis labels
mtext("number of contacted fragments", side=1, line=2.75, cex=0.95)
mtext("% baits", side=2, line=2.5, cex=0.95)

mtext(colnames(pc.degree.bait), at=apply(b, 2, mean), line=0.25, cex=0.75, side=1, las=2)

## legend & plot label
legend("topleft", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],  fill=dataset.colors[c("Original", "Simulated")], bty='n', cex=1, inset=c(0.05, -0.15), xpd=NA)
  
mtext("a", side=3, line=1, at=-8, font=2, cex=1.2)

##########################################################################

## bait degree

par(mar=c(4.5, 4.5, 2.1, 1.1))

ylim=c(0, max(as.numeric(pc.degree.frag)))

b=barplot(pc.degree.frag, beside=T, xlab='', ylim=ylim, space=c(0.4,1), names=rep("", dim(pc.degree.frag)[2]), ylab="", border=dataset.colors[c("Original", "Simulated")], col=dataset.colors[c("Original", "Simulated")], lwd=1.5,  mgp=c(3, 0.75, 0), cex.axis=1, las=2)

## axis labels
mtext("number of contacting baits", side=1, line=1.75, cex=0.95)
mtext("% fragments", side=2, line=2.5, cex=0.95)

mtext(colnames(pc.degree.frag), at=apply(b, 2, mean), line=0.25, cex=0.8, side=1, las=1)

## plot label
  
mtext("b", side=3, line=1, at=-8, font=2, cex=1.2)

##########################################################################

## correlation between number of observed contacts and number of genes in 500 kb window
## observed data

rho=round(cor(as.numeric(frag.degree.obs), stats.obs[names(frag.degree.obs), "nb_genes_500kb"], method="spearman"), digits=2)

par(mar=c(3.75, 4.5, 3, 1.1))
smoothScatter(as.numeric(frag.degree.obs), stats.obs[names(frag.degree.obs), "nb_genes_500kb"], pch=20, xlab="", ylab="", axes=F, cex=0.5, nbin=75)

axis(side=1, cex=1, mgp=c(3, 0.5, 0))
axis(side=2, cex=1, mgp=c(3, 0.75, 0), las=2)

mtext(paste("observed data, rho=",rho, sep=""), side=3, cex=0.85, line=0.25)

mtext("number of contacting baits", side=1, line=2, cex=0.95)
mtext("number of genes in 500 kb window", side=2, line=2.5, cex=0.95)

## plot label
  
mtext("c", side=3, line=1.5, at=-12, font=2, cex=1.2)

##########################################################################

## correlation between number of contacts and number of genes in 500 kb window
## simulated data

rho=round(cor(as.numeric(frag.degree.sim), stats.sim[names(frag.degree.sim), "nb_genes_500kb"], method="spearman"), digits=2)

par(mar=c(3.75, 4.5, 3, 1.1))
smoothScatter(as.numeric(frag.degree.sim), stats.sim[names(frag.degree.sim), "nb_genes_500kb"], pch=20, xlab="", ylab="", axes=F, cex=0.5, nbin=75)

axis(side=1, cex=1, mgp=c(3, 0.5, 0))
axis(side=2, cex=1, mgp=c(3, 0.75, 0), las=2)

mtext(paste("simulated data, rho=",rho, sep=""), side=3, cex=0.85, line=0.25)

mtext("number of contacting baits", side=1, line=2, cex=0.95)
mtext("number of genes in 500 kb window", side=2, line=2.5, cex=0.95)

## plot label
  
mtext("d", side=3, line=1.5, at=-6, font=2, cex=1.2)

##########################################################################


dev.off()

##########################################################################
