##############################################################################
setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")


source("parameters.R")

class_leg <- c("0",  "0.5",  "1", "1.5", "2")
par(lwd = 0.7)

pdf(paste(pathFigures, "SupplementaryFigure25.pdf", sep=""), width=7, height=3.5)

par(mai = c(0.65, 0.8, 0.5, 0.2)) # bottom, left, top, right
par(mfrow=c(1,2))

############### Supplementary Fig 25.A-B - Contact conservation by distance from TSS ##############

for (ref_sp in c("human", "mouse")){
  
  enhancers=enhancer.datasets[[ref_sp]]
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.contact.conservation.", ref_sp, ".Rdata", sep=""))
  
  if (ref_sp == "human"){YLIM=c(0,50); tg_sp="mouse"; lab="a"}else{YLIM=c(0, 50); label.enhancers=enhancers; tg_sp="human"; lab="b"}
  
  xlim=c(0.5, length(cons.dist[["ENCODE"]]["obs",])+0.5)
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM, xaxs="i", yaxs="i")
  
  for (enh in enhancer.datasets[[ref_sp]]){
    points(cons.dist[[enh]]["obs",],pch=20, col=col.enhancers[enh], cex=0.7)
    
    # Confidence intervals
    segments(x0=1:length(cons.dist.conf.low[[enh]]),y0=cons.dist.conf.low[[enh]]["obs",],
             x1=1:length(cons.dist.conf.low[[enh]]),y1=cons.dist.conf.high[[enh]]["obs",],
             col=col.enhancers[enh], lwd=0.3)
  }
  
  ## axis, legend & plot label
  axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0))
  mtext("Distance from TSS (Mb)", side=1, line=2)
  
  axis(side=2, mgp=c(3, 0.75, 0), las=2)
  mtext("% of conserved contacts", side=2, line=2.5)
  
  mtext(paste(ref_sp, "vs.", tg_sp), side=3, line=0.5)
  
  legend("bottomleft", col=col.enhancers, legend = label.enhancers, bty='n',pch=20, cex=0.8)
  mtext(lab, side=3, line=1, at=0.1, font=2)
}

############################################# 
dev.off()
