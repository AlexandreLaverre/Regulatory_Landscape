##############################################################################

source("../main_figures/parameters.R")

##############################################################################

pdf(paste(pathFigures, "SupplementaryMaterialFigure32.pdf", sep=""), width=6.85, height=3)
par(mfrow=c(1,2))
par(mar=c(3, 4, 2, 1)) # bottom, left, top, right

nb = 1
mtext.CEX=0.8

for (ref in c("human", "mouse")){
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.", ref, ".stats.RData", sep=""))
  tg=setdiff(c("human", "mouse"), ref)
  
  ## conservation by distance
  cons.dist = cons.dist[["all"]]
  cons.dist.conf.low = cons.dist.conf.low[["all"]]
  cons.dist.conf.high = cons.dist.conf.high[["all"]]
  
  YLIM=c(-1,50)
  class_leg <- c("0.05",  "0.5",  "1", "1.5", "2")
  
  xlim=c(0.5, length(cons.dist[["ENCODE"]]["obs",])+0.5)
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM,  yaxs="i")

  for (enh in enhancer.datasets[[ref]]){
    nbclass=length(cons.dist.conf.low[[enh]]["obs",])
    
    points(cons.dist[[enh]]["obs",],pch=21, col=col.enhancers[enh], bg=col.enhancers[enh], cex=0.4)
    segments(1:nbclass, cons.dist.conf.low[[enh]]["obs",], 1:nbclass, cons.dist.conf.high[[enh]]["obs",], col=col.enhancers[enh], cex=0.5)
    
    points(cons.dist[[enh]]["sim",],pch=21, col=col.enhancers[enh], bg="white")
    segments(1:nbclass, cons.dist.conf.low[[enh]]["sim",], 1:nbclass, cons.dist.conf.high[[enh]]["sim",], col=col.enhancers[enh])
  }
 
  ## axis, legend & plot label
  axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=0.8)
  mtext("distance from promoter region (Mb)", side=1, line=2, cex=mtext.CEX)
  
  axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=0.8)
  mtext("% conserved contacts", side=2, line=2.5,  cex=mtext.CEX)
  
  mtext(paste(ref, 'vs.', tg, sep=" "), side=3, line=0.5, at=22, cex=0.8)

  ## plot label
  mtext(letters[nb], side=3, line=1, at=-11, font=2, cex=1.1)
  
  par(xpd=TRUE)
  
  if (nb == 1){
    legend("topleft", col=col.enhancers, legend = label.enhancers, box.col="white", bg="white", pch=20, cex=0.75, inset=c(0.01, -0.05))
  } else{
     legend("topright", legend=c("PCHi-C data", "simulated data"),  pch=21, pt.bg=c("black", "white"),  bty='n', inset=c(-0.05, 0), xpd=NA, cex=0.75)
  }
    
  nb = nb+1
    
}

dev.off()

##############################################################################
