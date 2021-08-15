###############################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

##############################################################################

if(load){
 
  ## expression divergence
  
  load(paste(pathFigures, "RData/data.human.CM2019.AllOrgans.expdiv.RData", sep=""))
  expdiv.human=expdiv

  load(paste(pathFigures, "RData/data.mouse.CM2019.AllOrgans.expdiv.RData", sep=""))
  expdiv.mouse=expdiv
  
  rm("expdiv")

  ## number of enhancers
  
  enh="ENCODE"

  load(paste(pathFigures, "RData/data.neighbor.enhancers.RData", sep="")) ## for now we keep all distance ranges

  neighbors.human=gene.enhancer.contacts[["human"]][[enh]][["real"]]
  neighbors.mouse=gene.enhancer.contacts[["mouse"]][[enh]][["real"]]

  ## predicted regulatory regions
  regions.human=gene.enhancer.contacts[["human"]][["regions"]]
  regions.mouse=gene.enhancer.contacts[["mouse"]][["regions"]]

  load=FALSE
}

##############################################################################

if(prepare){
 
  ## number of enhancers

  nb.enh.human=as.numeric(table(as.factor(neighbors.human$GeneID)))
  names(nb.enh.human)=levels(as.factor(neighbors.human$GeneID))
  
  nb.enh.mouse=as.numeric(table(neighbors.mouse$GeneID))
  names(nb.enh.mouse)=levels(as.factor(neighbors.mouse$GeneID))

  ## same order
  nb.enh.human=nb.enh.human[rownames(expdiv.human)]
  nb.enh.mouse=nb.enh.mouse[rownames(expdiv.mouse)]

  ## regions
  rownames(regions.human)=regions.human$gene_id
  rownames(regions.mouse)=regions.mouse$gene_id

  regions.human$size=regions.human$end_region-regions.human$start_region+1
  regions.mouse$size=regions.mouse$end_region-regions.mouse$start_region+1

  regions.human=regions.human[rownames(expdiv.human),]
  regions.mouse=regions.mouse[rownames(expdiv.mouse),]

  ## enhancer number class

  class.enh.human=cut(nb.enh.human, breaks=c(1, 10, 20, 30, 40, max(nb.enh.human, na.rm=T)), include.lowest=T)
  names(class.enh.human)=names(nb.enh.human)
    
  class.enh.mouse=cut(nb.enh.mouse, breaks=c(1, 10, 20, 30, 40, max(nb.enh.mouse, na.rm=T)), include.lowest=T)
  names(class.enh.mouse)=names(nb.enh.mouse)

  enh.classes=c("1-10", "11-20", "21-30", "31-40", ">40")

  ## region size

  size.breaks=c(0, 50e3, 100e3, 150e3, 200e3, 2e6+1)  
  class.size.human=cut(regions.human$size, breaks=size.breaks, include.lowest=T)
  names(class.size.human)=rownames(regions.human)
  
  class.size.mouse=cut(regions.mouse$size, breaks=size.breaks, include.lowest=T)
  names(class.size.mouse)=rownames(regions.mouse)

  size.classes=c("0-50", "50-100", "100-150", "150-200", ">200")
  
  prepare=FALSE
}

##############################################################################

prepare.plot <-function(y, class){
  lev=levels(class)

  medians=c()
  ci.low=c()
  ci.high=c()

  for(l in lev){
    g=which(class==l)
    medians=c(medians, median(y[g], na.rm=T))
    b=boxplot(y[g], plot=F)
    ci.low=c(ci.low, b$conf[1,1])
    ci.high=c(ci.high, b$conf[2,1])
  }

  names(medians)=lev
  names(ci.low)=lev
  names(ci.high)=lev
  
  ylim=range(c(ci.low, ci.high))
  diffy=diff(ylim)/20
  ylim=ylim+c(-diffy, diffy)
  
  xlim=c(1, length(lev))

  results=list("xlim"=xlim, "ylim"=ylim, "medians"=medians, "ci.low"=ci.low, "ci.high"=ci.high)

  return(results)
}

##############################################################################

pdf(file=paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_neighbor_enhancers_expression_patterns.pdf", sep=""), width=7, height=6.25)

##############################################################################

m=matrix(rep(NA, 8), nrow=2)

m[1,]=c(rep(1, 2), rep(2,2))
m[2,]=c(rep(3, 2), rep(4,2))

layout(m)

tinyx=0.15
cex.axis=1.1
cex.mtext=0.85
cex.label=1.1

#############################################################################################

## first plot: mean expression level against number of contacts

plot.human=prepare.plot(expdiv.human$MeanRPKM, class.enh.human)
plot.mouse=prepare.plot(expdiv.mouse$MeanRPKM, class.enh.mouse)

ylim=range(c(plot.human[["ylim"]], plot.mouse[["ylim"]]))
xlim=plot.human[["xlim"]]+c(-0.5, 0.5)

xpos=1:length(enh.classes)

par(mar=c(4.1, 4.5, 2.75, 0.75))

plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
points(xpos-tinyx, plot.human[["medians"]], pch=20, col="darkred")
segments(xpos-tinyx, plot.human[["ci.low"]], xpos-tinyx, plot.human[["ci.high"]], col="darkred")

points(xpos+tinyx, plot.mouse[["medians"]], pch=20, col="navy")
segments(xpos+tinyx, plot.mouse[["ci.low"]], xpos+tinyx, plot.mouse[["ci.high"]], col="navy")

axis(side=1, at=xpos, labels=enh.classes, cex.axis=cex.axis, mgp=c(3, 0.75, 0))
axis(side=2, cex.axis=cex.axis, mgp=c(3, 0.5, 0))

mtext("mean expression level (RPKM)", side=2, line=2.25, cex=cex.mtext)
mtext("nb. neighbor enhancers", side=1, line=2.25, cex=cex.mtext)

mtext("A", side=3, font=2, cex=cex.label, at=-0.45, line=1)

#############################################################################################

## second plot: mean expression specificity against number of contacts

plot.human=prepare.plot(expdiv.human$MeanTau, class.enh.human)
plot.mouse=prepare.plot(expdiv.mouse$MeanTau, class.enh.mouse)

ylim=range(c(plot.human[["ylim"]], plot.mouse[["ylim"]]))
xlim=plot.human[["xlim"]]+c(-0.5, 0.5)

xpos=1:length(enh.classes)

par(mar=c(4.1, 4.5, 2.75, 0.75))

plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
points(xpos-tinyx, plot.human[["medians"]], pch=20, col="darkred")
segments(xpos-tinyx, plot.human[["ci.low"]], xpos-tinyx, plot.human[["ci.high"]], col="darkred")

points(xpos+tinyx, plot.mouse[["medians"]], pch=20, col="navy")
segments(xpos+tinyx, plot.mouse[["ci.low"]], xpos+tinyx, plot.mouse[["ci.high"]], col="navy")

axis(side=1, at=xpos, labels=enh.classes, cex.axis=cex.axis, mgp=c(3, 0.75, 0))
axis(side=2, cex.axis=cex.axis, mgp=c(3, 0.5, 0))

mtext("mean expression specificity", side=2, line=2.25, cex=cex.mtext)
mtext("nb. neighbor enhancers", side=1, line=2.25, cex=cex.mtext)

mtext("B", side=3, font=2, cex=cex.label, at=-0.45, line=1)

#############################################################################################

## third plot: expression conservation against number of contacts

plot.human=prepare.plot(expdiv.human$CorrectedSpearman, class.enh.human)
plot.mouse=prepare.plot(expdiv.mouse$CorrectedSpearman, class.enh.mouse)

ylim=range(c(plot.human[["ylim"]], plot.mouse[["ylim"]]))
xlim=plot.human[["xlim"]]+c(-0.5, 0.5)

xpos=1:length(enh.classes)

par(mar=c(4.1, 4.5, 2.75, 0.75))

plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
points(xpos-tinyx, plot.human[["medians"]], pch=20, col="darkred")
segments(xpos-tinyx, plot.human[["ci.low"]], xpos-tinyx, plot.human[["ci.high"]], col="darkred")

points(xpos+tinyx, plot.mouse[["medians"]], pch=20, col="navy")
segments(xpos+tinyx, plot.mouse[["ci.low"]], xpos+tinyx, plot.mouse[["ci.high"]], col="navy")

axis(side=1, at=xpos, labels=enh.classes, cex.axis=cex.axis, mgp=c(3, 0.75, 0))
axis(side=2, cex.axis=cex.axis, mgp=c(3, 0.5, 0))

mtext("expression conservation", side=2, line=2.25, cex=cex.mtext)
mtext("nb. neighbor enhancers", side=1, line=2.25, cex=cex.mtext)

mtext("C", side=3, font=2, cex=cex.label, at=-0.45, line=1)

#############################################################################################

## fourth plot: expression conservation against region size

plot.human=prepare.plot(expdiv.human$CorrectedSpearman, class.size.human)
plot.mouse=prepare.plot(expdiv.mouse$CorrectedSpearman, class.size.mouse)

ylim=range(c(plot.human[["ylim"]], plot.mouse[["ylim"]]))
xlim=plot.human[["xlim"]]+c(-0.5, 0.5)

xpos=1:length(enh.classes)

par(mar=c(4.1, 4.5, 2.75, 0.75))

plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
points(xpos-tinyx, plot.human[["medians"]], pch=20, col="darkred")
segments(xpos-tinyx, plot.human[["ci.low"]], xpos-tinyx, plot.human[["ci.high"]], col="darkred")

points(xpos+tinyx, plot.mouse[["medians"]], pch=20, col="navy")
segments(xpos+tinyx, plot.mouse[["ci.low"]], xpos+tinyx, plot.mouse[["ci.high"]], col="navy")

axis(side=1, at=xpos[c(1,2, 4, 5)], labels=size.classes[c(1,2, 4, 5)], cex.axis=cex.axis, mgp=c(3, 0.75, 0))
axis(side=1, at=xpos[3], labels=size.classes[3], cex.axis=cex.axis, mgp=c(3, 0.75, 0))
axis(side=2, cex.axis=cex.axis, mgp=c(3, 0.5, 0))

mtext("expression conservation", side=2, line=2.25, cex=cex.mtext)
mtext("distance to next TSS (kb)", side=1, line=2.25, cex=cex.mtext)

mtext("D", side=3, font=2, cex=cex.label, at=-0.45, line=1)

#############################################################################################

dev.off()

#############################################################################################
