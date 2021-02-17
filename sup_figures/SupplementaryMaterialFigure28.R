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
  
  class_leg <- c("0.05",  "0.5",  "1", "1.5", "2")
   
  enh="ENCODE"
  
  load=FALSE
}

##############################################################################

pdf(paste(pathFigures, "SupplementaryMaterialFigure25.pdf", sep=""), width=6.85, height=3.5)

par(mar=c(3.5, 3.25, 2.1, 1.1))
par(mfrow=c(1,2))
par(lwd = 0.7)

############################################################################## 

for (ref_sp in c("human", "mouse")){
  
  tg=setdiff(c("human", "mouse"), ref_sp)
  
  enhancers=enhancer.datasets[[ref_sp]]
  label.enhancers=enhancers
  
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.", ref_sp, ".stats.RData", sep=""))
  
  if (ref_sp == "human"){
    YLIM=c(-1,80)
    lab="a"
  }else{
    YLIM=c(-1, 80)
    lab="b"
  }
  
  xlim=c(-0.5, length(cons.dist[["all"]][[enh]]["obs",])+1)
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM, xaxs="i", yaxs="i")
  
  for (data in c("obs", "sim")){
    
    if (data == "obs"){
      data.name = "Original"
    }else{
      data.name = "Simulated"
    }
    
    for (onto in ontologies){
      if (data == "obs"){
        points(cons.dist[[onto]][[enh]][data,], pch=21, bg=col.ontologies[onto], col=col.ontologies[onto], cex=0.7)
      } else{
        points(cons.dist[[onto]][[enh]][data,], pch=21, bg="white",  col=col.ontologies[onto], cex=0.7)
      }
      
      # Confidence intervals
      segments(x0=1:length(cons.dist.conf.low[[onto]][[enh]]),y0=cons.dist.conf.low[[onto]][[enh]][data,],
               x1=1:length(cons.dist.conf.low[[onto]][[enh]]),y1=cons.dist.conf.high[[onto]][[enh]][data,],
               col=col.ontologies[onto], lwd=0.3)
    }
  }
  
  ## axis, legend & plot label
  axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=0.8)
  mtext("distance from promoter region (Mb)", side=1, line=2, cex=0.8)
  
  axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=0.8)
  mtext("% conserved contacts", side=2, line=2.25, cex=0.8)
  
  mtext(paste(ref_sp, "vs.", tg), side=3, line=0.5, cex=0.8)
  
  if (ref_sp == "human"){
    legend("topleft", col=col.ontologies, legend = c("developmental genes", "immune genes", "other genes"), bty='n', pch=20, cex=0.75)

  } else{
    legend("topleft", legend=c("PCHi-C data", "simulated data"),  pch=c(21,21), pt.bg=c("black", "white"),
           bty='n', xpd=NA, cex=0.75)
  }
  
  mtext(lab, side=3, line=1, at=-9, font=2, cex=1.1)

}

dev.off()

############################################################################## 
