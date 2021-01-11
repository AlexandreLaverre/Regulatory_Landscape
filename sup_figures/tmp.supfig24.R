##############################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  
  source("../main_figures/parameters.R")
}

##############################################################################

if(load){
  ref = "mouse"
  tg=setdiff(c("human", "mouse"), ref)
  
  enhancers=enhancer.datasets[[ref]]
  
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData", sep=""))
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.", ref, ".stats.RData", sep=""))

  enh="ENCODE"
    
  load=FALSE
}

##############################################################################

if(prepare){

  samples.ref=sampleinfo[[ref]][,"Sample.ID"]
  samples.tg=sampleinfo[[tg]][,"Sample.ID"]

  cells.ref=sampleinfo[[ref]][,"Broad.cell.type.or.tissue"]
  cells.tg=sampleinfo[[tg]][,"Broad.cell.type.or.tissue"]

  names(cells.ref)=samples.ref
  names(cells.tg)=samples.tg

  nbmax=length(unique(cells.tg))

  ## observed and simulated contacts
    
  cc.obs=contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]][["obsobs"]]
  cc.sim=contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]][["simsim"]]

  cc.obs$target_nb_cell=apply(cc.obs[,samples.tg],1, function(x) length(unique(cells.tg[which(x>0)])))
  cc.sim$target_nb_cell=apply(cc.sim[,samples.tg],1, function(x) length(unique(cells.tg[which(x>0)])))
  
  cc.obs$ref_nb_cell=apply(cc.obs[,samples.ref],1, function(x) length(unique(cells.ref[which(x>0)])))
  cc.sim$ref_nb_cell=apply(cc.sim[,samples.ref],1, function(x) length(unique(cells.ref[which(x>0)])))

  ## conservation by distance class

  cc.obs$dist_class=cut(cc.obs$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest=T)
  cc.sim$dist_class=cut(cc.sim$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest=T)

  prop.obs.dist=100*tapply(cc.obs$target_nb_cell, cc.obs$dist_class, function(x) length(which(x>0))/length(x))
  prop.sim.dist=100*tapply(cc.sim$target_nb_cell, cc.sim$dist_class, function(x) length(which(x>0))/length(x))

  ci.low.obs.dist=100*tapply(cc.obs$target_nb_cell, cc.obs$dist_class, function(x) prop.test(length(which(x>0)), length(x))$conf.int[1])
  ci.low.sim.dist=100*tapply(cc.sim$target_nb_cell, cc.sim$dist_class, function(x) prop.test(length(which(x>0)), length(x))$conf.int[1])

  ci.high.obs.dist=100*tapply(cc.obs$target_nb_cell, cc.obs$dist_class, function(x) prop.test(length(which(x>0)), length(x))$conf.int[2])
  ci.high.sim.dist=100*tapply(cc.sim$target_nb_cell, cc.sim$dist_class, function(x) prop.test(length(which(x>0)), length(x))$conf.int[2])

  ## conservation by number of cell types

  nbmax=max(cc.obs$ref_nb_cell)
  cc.obs$cell_class=cut(cc.obs$ref_nb_cell, breaks=c(seq(from=0, to=5, by=1), nbmax), include.lowest=T)
  cc.sim$cell_class=cut(cc.sim$ref_nb_cell, breaks=c(seq(from=0, to=5, by=1), nbmax), include.lowest=T)

  prop.obs.cell=100*tapply(cc.obs$target_nb_cell, cc.obs$cell_class, function(x) length(which(x>0))/length(x))
  prop.sim.cell=100*tapply(cc.sim$target_nb_cell, cc.sim$cell_class, function(x) length(which(x>0))/length(x))

  ci.low.obs.cell=100*tapply(cc.obs$target_nb_cell, cc.obs$cell_class, function(x) prop.test(length(which(x>0)), length(x))$conf.int[1])
  ci.low.sim.cell=100*tapply(cc.sim$target_nb_cell, cc.sim$cell_class, function(x) prop.test(length(which(x>0)), length(x))$conf.int[1])

  ci.high.obs.cell=100*tapply(cc.obs$target_nb_cell, cc.obs$cell_class, function(x) prop.test(length(which(x>0)), length(x))$conf.int[2])
  ci.high.sim.cell=100*tapply(cc.sim$target_nb_cell, cc.sim$cell_class, function(x) prop.test(length(which(x>0)), length(x))$conf.int[2])

  matrix.cell=matrix(c(prop.obs.cell, prop.sim.cell), nrow=2, byrow=T)

  ## conservation in comparable cell types
  
  prepare=FALSE
}

##############################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##############################################################################

pdf(paste(pathFigures, "SupplementaryFigure24.pdf", sep=""), width=6.85, height=2.5)

m=matrix(rep(NA, 1*13), nrow=1)
m[1,]=c(rep(1,4), rep(2,3), rep(3,6))
layout(m)

mtext.CEX=0.75
par(lwd = 1.5)

#########################################################################################

## conservation as a function of nb cell types

par(mar=c(3.1, 4.5, 2.5, 0.25))
b=barplot(matrix.cell, border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")], mgp=c(3, 0.75, 0), las=2, beside=T, cex.axis=1.05, ylim=c(0, 70), space=c(0.2, 1))

segments(b[1,], ci.low.obs.cell, b[1,], ci.high.obs.cell)
segments(b[2,], ci.low.sim.cell, b[2,], ci.high.sim.cell)

axis(side=1, at=apply(b,2,mean), labels=c(as.character(1:5), "6+"), cex.axis=1.05, mgp=c(3, 0.75, 0))
mtext("% conserved contacts", side=2, line=3, cex=mtext.CEX)
mtext("nb. cell types", side=1, line=2.0, cex=mtext.CEX)

mtext("a", side=3, line=1.25, at=-5.65, font=2, cex=1.05)

## legend 
legend("topleft", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],
       fill=dataset.colors[c("Original", "Simulated")], bty='n', 
       inset=c(0.05, -0.12), xpd=NA)

####################################################################################

## conservation in common cell types

YMAX=25

par(mar=c(4.1, 4.5, 2.5, 1.5))

b=barplot(cons.common.cell[[enh]], beside=T, names=rep("", dim(cons.common.cell[[enh]])[2]), ylim=c(0,YMAX), space=c(0.2,1),
          border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
          mgp=c(3, 0.75, 0), cex.axis=1.05, las=2)

arrows(x0=b,y0=cons.common.cell.conf.low[[enh]],y1=cons.common.cell.conf.high[[enh]],angle=90,code=3,length=0.01)

## axis labels

mtext("% conserved contacts", side=2, line=2.5, cex=mtext.CEX)

label.cells = c("ESC", "adipo", "B cells")
mtext(label.cells, at=apply(b, 2, mean), side=1, line=0.5, cex=mtext.CEX, las=2)


mtext("b", side=3, line=1, at=-4, font=2, cex=1.05)

#########################################################################################

## conservation as a function of distance

xpos=1:length(prop.obs.dist)
xlim=c(-1, max(xpos)+1)
ylim=c(-0.5, max(c(ci.high.obs.dist, ci.high.sim.dist))+5)

par(mar=c(4.1, 4.5, 2.5, 1.5))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

points(xpos, prop.obs.dist, col=dataset.colors["Original"], pch=20)
segments(xpos, ci.low.obs.dist, xpos, ci.high.obs.dist, col=dataset.colors["Original"])

points(xpos, prop.sim.dist, col=dataset.colors["Simulated"], pch=20)
segments(xpos, ci.low.sim.dist, xpos, ci.high.sim.dist, col=dataset.colors["Simulated"])

axis(side=1, at=c(xpos[1], xpos[10], xpos[20], xpos[30], xpos[39]+1), labels=c(0.025, 0.5, 1, 1.5, 2), cex.axis=1.05)
axis(side=2, cex.axis=1.05, las=2)
mtext("% conserved contacts", side=2, line=3, cex=mtext.CEX)
mtext("distance to promoters (Mb)", side=1, line=2.5, cex=mtext.CEX)

mtext("c", side=3, line=1.25, at=-9.1, font=2, cex=1.05)

#########################################################################################

dev.off()

####################################################################################
