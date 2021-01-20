##############################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  
  source("../main_figures/parameters.R")
}

##############################################################################
if(load){
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData", sep=""))
  
  ontologies = c("dvpt", "immune", "other")
  col.ontologies = c("forestgreen", "orange", "red")
  names(col.ontologies) = ontologies
  
  class_leg <- c("0",  "0.5",  "1", "1.5", "2")
  par(lwd = 0.7)
  
  enh="ENCODE"
  load=FALSE
}


##############################################################################

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

##############################################################################

#pdf(paste(pathFigures, "SupplementaryFigure24.pdf", sep=""), width=6.85, height=6)

par(mai = c(0.65, 0.8, 0.5, 0.2)) # bottom, left, top, right
par(mfrow=c(1,2))

############################################################################## 

for (ref_sp in c("human", "mouse")){
  
  tg=setdiff(c("human", "mouse"), ref_sp)
  
  enhancers=enhancer.datasets[[ref_sp]]
  label.enhancers=enhancers
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.", ref_sp, ".stats.Rdata", sep=""))
  
  if (ref_sp == "human"){YLIM=c(-1,80); lab="a"}else{YLIM=c(-1, 80); lab="b"}

  xlim=c(0.5, length(cons.dist[["all"]][[enh]]["obs",])+0.5)
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM, xaxs="i", yaxs="i")
  
  for (data in c("obs", "sim")){
    if (data == "obs"){data.name = "Original"}else{data.name = "Simulated"}
    
    for (onto in ontologies){
      if (data == "obs"){points(cons.dist[[onto]][[enh]][data,], type="l", col=col.ontologies[onto], cex=0.7)
        }else{points(cons.dist[[onto]][[enh]][data,], pch=20, col=col.ontologies[onto], cex=0.7)}
      
      # Confidence intervals
      segments(x0=1:length(cons.dist.conf.low[[onto]][[enh]]),y0=cons.dist.conf.low[[onto]][[enh]][data,],
               x1=1:length(cons.dist.conf.low[[onto]][[enh]]),y1=cons.dist.conf.high[[onto]][[enh]][data,],
               col=col.ontologies[onto], lwd=0.3)
    }
  }
  
  ## axis, legend & plot label
  axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0))
  mtext("Distance from TSS (Mb)", side=1, line=2)
  
  axis(side=2, mgp=c(3, 0.75, 0), las=2)
  mtext("% of conserved contacts", side=2, line=2.5)
  
  mtext(paste(ref_sp, "vs.", tg), side=3, line=0.5)
  
  legend("topleft", col=col.ontologies, legend = c("developmental", "immune", "other"), bty='n', lty=1)
  mtext(lab, side=3, line=1, at=0.1, font=2)

}

dev.off()