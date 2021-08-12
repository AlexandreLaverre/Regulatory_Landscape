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
  sp="human"
  enh="ENCODE"

  print("loading data")

  load(paste(pathFigures, "RData/data.neighbor.enhancers.RData", sep=""))
  neighbors=gene.enhancer.contacts[[sp]][[enh]][["real"]]

  load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
  contacts=gene.enhancer.contacts[[sp]][[enh]][["real"]]

  print("done")
  
  load=FALSE
}

#############################################################################################

if(prepare){
  
  ## select common genes

  common.genes=intersect(contacts$gene, neighbors$GeneID)

  contacts.common=contacts[which(contacts$gene%in%common.genes),]
  neighbors.common=neighbors[which(neighbors$GeneID%in%common.genes),]

  ## number of enhancers per gene

  nb.contacts.common=as.numeric(table(as.factor(contacts$gene)))
  names(nb.contacts.common)=levels(as.factor(contacts$gene))
  
  nb.neighbors.common=as.numeric(table(as.factor(neighbors$GeneID)))
  names(nb.neighbors.common)=levels(as.factor(neighbors$GeneID))

  nb.neighbors.common=nb.neighbors.common[common.genes]
  nb.contacts.common=nb.contacts.common[common.genes]

  ## shared enhancers

  nb.shared.common=unlist(lapply(common.genes, function(x) length(intersect(neighbors$EnhancerID[which(neighbors$GeneID==x)], contacts$enhancer[which(contacts$gene==x)]))))
  names(nb.shared.common)=common.genes

  fr.common.neighbors=100*nb.shared.common/nb.neighbors.common
  fr.common.contacts=100*nb.shared.common/nb.contacts.common
  

  ## median distance per gene

  median.dist.neighbors=tapply(neighbors.common$Distance, as.factor(neighbors.common$GeneID), median)
  median.dist.contacts=tapply(contacts.common$dist, as.factor(contacts.common$gene), median)
  median.dist.contacts=median.dist.contacts[common.genes]
  median.dist.neighbors=median.dist.neighbors[common.genes]
  
  prepare=FALSE
}

#############################################################################################

pdf(file=paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_neighbor_enhancers.pdf", sep=""), width=6.85, height=3.1)

m=matrix(rep(NA, 10), nrow=1)

m[1,]=c(rep(1, 4), rep(2,2), rep(3, 4))

layout(m)

#############################################################################################

## number of enhancers, neighbor vs. contacted

lim=c(0, max(c(nb.contacts.common, nb.neighbors.common)))

par(mar=c(4.5, 3.5, 2.5, 1.1))

smoothScatter(nb.contacts.common, nb.neighbors.common, pch=20, xlab="", ylab="", axes=F, xlim=lim, ylim=lim)

abline(0, 1, col="black")

mtext("nb. contacted enhancers", side=1, line=2.1, cex=0.85)
mtext("nb. neighbor enhancers", side=2, line=2.1, cex=0.85)

rho=cor(nb.contacts.common, nb.neighbors.common, method="spearman")

mtext(paste("rho =",round(rho, digits=2)), side=3, cex=0.85, line=0.25) 

axis(side=1, cex.axis=0.95, mgp=c(3, 0.5, 0))
axis(side=2, cex.axis=0.95, mgp=c(3, 0.5, 0))

higher=100*length(which(nb.contacts.common>nb.neighbors.common))/length(nb.contacts.common)
lower=100*length(which(nb.contacts.common<=nb.neighbors.common))/length(nb.contacts.common)

text(paste(round(higher, digits=1), "%", sep=""), x=330, y=10, cex=1.1)
text(paste(round(lower, digits=1), "%", sep=""), x=25, y=350, cex=1.1)

mtext("A", side=3, at=-75, cex=1.1, font=2, line=1.1)

#############################################################################################

par(mar=c(5.5, 3.5, 2.5, 1.1))
boxplot(nb.contacts.common, nb.neighbors.common, col="white", border=c("darkorange", "black"), boxwex=0.5, axes=F, outline=F)

axis(side=1, cex.axis=0.95, mgp=c(3, 0.5, 0), at=c(1,2), labels=rep("", 2))
axis(side=2, cex.axis=0.95, mgp=c(3, 0.5, 0))

mtext("nb. enhancers", side=2, line=2.1, cex=0.9)
mtext(c("contacts", "neighbors"), at=c(1, 2), side=1, cex=0.85, las=2, line=0.55)

mtext("B", side=3, at=-0.5, cex=1.1, font=2, line=1.1)

#############################################################################################

## distance distribution, all gene-enhancer pairs

d.n=density(neighbors.common$Distance/1e3)
d.c=density(contacts.common$dist/1e3, bw=d.n$bw)

xlim=range(c(d.n$x, d.c$x))
ylim=range(c(d.n$y, d.c$y))

par(mar=c(4.5, 3.5, 2.5, 1.1))

plot(d.n$x, d.n$y, type="l", col="black", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F, lwd=1.25)
lines(d.c$x, d.c$y, type="l", col="darkorange", lwd=1.25)

mtext("gene-enhancer distance (kb)", side=1, line=2.1, cex=0.85)
mtext("density", side=2, line=2.1, cex=0.85)

axis(side=1, cex.axis=0.95, mgp=c(3, 0.5, 0))
axis(side=2, cex.axis=0.95, mgp=c(3, 0.5, 0))

median.neighbors=median(neighbors.common$Distance, na.rm=T)/1e3
median.contacts=median(contacts.common$dist, na.rm=T)/1e3

abline(v=median.neighbors, col="black", lty=2)
abline(v=median.contacts, col="darkorange", lty=2)

mtext(round(median.neighbors, digits=0), side=3, at=median.neighbors, las=1, cex=0.75)
mtext(round(median.contacts, digits=0), side=3, at=median.contacts, las=1, cex=0.75, col="darkorange")

mtext("median", line=1, cex=0.75, col="black", at=mean(c(median.neighbors, median.contacts)), side=3)

legend("topright", inset=-0.05, col=c("black", "darkorange"), lty=1, legend=c("neighbor enhancers", "contacted enhancers"), cex=1.1, bty="n", xpd=NA)

mtext("C", side=3, at=-480, cex=1.1, font=2, line=1.1)

#############################################################################################

dev.off()

#############################################################################################

