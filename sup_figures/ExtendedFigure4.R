#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  
  source("../main_figures/parameters.R")

  ref="human"
  tg="mouse"
  
  enh="ENCODE"
}

##########################################################################

if(load){

  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))

  stats.obs=enhancer.statistics[[ref]][[enh]][["original"]]
  stats.sim=enhancer.statistics[[ref]][[enh]][["simulated"]]

  stats.obs=stats.obs[which(stats.obs$all_exon_bp==0),]
  stats.obs=stats.obs[which(stats.obs$all_exon_bp==0),]

  load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref,"2", tg,".RData", sep=""))

  stats.obs$pcungapped=pcungapped[rownames(stats.obs)]
  stats.sim$pcungapped=pcungapped[rownames(stats.sim)]
  
  load=FALSE
}

##########################################################################

if(prepare){

  class.nbgenes.obs=cut(stats.obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(stats.obs$nb_genes_500kb)), include.lowest=T)
  class.nbgenes.sim=cut(stats.sim$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(stats.sim$nb_genes_500kb)), include.lowest=T)

  mean.cons.obs=tapply(stats.obs$pcungapped, class.nbgenes.obs, mean, na.rm=T)
  ci.low.obs=tapply(stats.obs$pcungapped, class.nbgenes.obs, function(x) t.test(x)[["conf.int"]][1])
  ci.high.obs=tapply(stats.obs$pcungapped, class.nbgenes.obs, function(x) t.test(x)[["conf.int"]][2])
  
  mean.cons.sim=tapply(stats.sim$pcungapped, class.nbgenes.sim, mean, na.rm=T)
  ci.low.sim=tapply(stats.sim$pcungapped, class.nbgenes.sim, function(x) t.test(x)[["conf.int"]][1])
  ci.high.sim=tapply(stats.sim$pcungapped, class.nbgenes.sim, function(x) t.test(x)[["conf.int"]][2])
    
  prepare=F
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "ExtendedFigure4b.pdf", sep=""), width=4.49, height=4)

m=matrix(c(rep(1,2), rep(2, 4)), nrow=1)
layout(m)

###########################################################################
## nb genes

par(mar=c(6.75, 4.1, 2.1, 1.1))
B(stats.obs$nb_genes_500kb, stats.sim$nb_genes_500kb, pch=20, col="white", border=dataset.colors, outline=F, boxwex=0.5, notch=T)
mtext("number of genes within 500kb", side=2, line=3, cex=0.75)

mtext(c("PCHi-C data", "simulated data"), at=1:2, side=1, line=0.5, cex=0.75, las=2)

###########################################################################

par(mar=c(6.75, 4.1, 2.1, 1.1))
xpos=1:length(mean.cons.obs)
xlim=range(xpos)+c(-0.5, 0.5)
ylim=range(c(ci.low.obs, ci.high.obs, ci.low.sim, ci.high.sim))

smallx=c(-0.1, 0.1)

plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim)

points(xpos+smallx[1], mean.cons.obs, col=dataset.colors["Original"], pch=20)
segments(xpos+smallx[1], ci.low.obs, xpos+smallx[1], ci.high.obs, col=dataset.colors[["Original"]])

points(xpos+smallx[2], mean.cons.sim, col=dataset.colors["Simulated"], pch=20)
segments(xpos+smallx[2], ci.low.sim, xpos+smallx[2], ci.high.sim, col=dataset.colors[["Simulated"]])

abline(v=xpos[-length(xpos)]+0.5, lty=3, col="gray40")

axis(side=2, mgp=c(3, 0.75, 0), las=2)
axis(side=1, mgp=c(3, 0.75, 0), at=xpos, labels=rep("", length(xpos)))

mtext(names(mean.cons.obs), at=xpos, side=1, line=1, las=2, cex=0.75)

mtext("number of genes within 500 kb", side=1, line=4.5, cex=0.75)
mtext("sequence conservation", side=2, line=3, cex=0.75)

legend("topright", col=dataset.colors, legend = c("PCHi-C data", "simulated data"), box.col="white", bg="white", pch=20, cex=1.1, inset=c(0.01,-0.01), xpd=T)

###########################################################################


dev.off()

###########################################################################

## que se passe-t-il sur les rearrangements ?
## Sequences au hasard au lieu des enhancers
## prendre que enhancers avec 0 chevauchement exonique

## activité des enhancers par rapport à l'activité du gène contacté
## synténie : comment calculer la densité en gènes ? 

## 
