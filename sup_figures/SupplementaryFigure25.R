##############################################################################
setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")


source("parameters.R")

class_leg <- c("0",  "0.5",  "1", "1.5", "2")
par(lwd = 0.7)

pdf(paste(pathFigures, "SupplementaryFigure25.pdf", sep=""), width=6.85, height=6)

par(mai = c(0.65, 0.8, 0.5, 0.2)) # bottom, left, top, right
par(mfrow=c(2,2))

############### Supplementary Fig 25.A-B - Contact conservation by distance from TSS ##############
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

############################################# 
ontologies = c("dvpt", "immune", "other")
col.ontologies = c("forestgreen", "orange", "red")
names(col.ontologies) = ontologies


for (ref_sp in c("mouse", "human")){
  
  enhancers=enhancer.datasets[[ref_sp]]
  label.enhancers=enhancers
  load(paste(pathFigures, "RData/data.contact.conservation.", ref_sp, ".Rdata", sep=""))
  
  if (ref_sp == "human"){YLIM=c(0,50); tg_sp="mouse"; lab="a"
  }else{YLIM=c(0, 50); tg_sp="human"; lab="b"}

  xlim=c(0.5, length(cons.dist[["all"]][["ENCODE"]]["obs",])+0.5)
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM, xaxs="i", yaxs="i")
  
  for (enh in enhancer.datasets[[ref_sp]]){
    points(cons.dist[["all"]][[enh]]["obs",],pch=20, col=col.enhancers[enh], cex=0.5)
    
    # Confidence intervals
    segments(x0=1:length(cons.dist.conf.low[["all"]][[enh]]),y0=cons.dist.conf.low[["all"]][[enh]]["obs",],
             x1=1:length(cons.dist.conf.low[["all"]][[enh]]),y1=cons.dist.conf.high[["all"]][[enh]]["obs",],
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


########################################################################################## 
YMAX=40

par(lwd = 1.5)

b=barplot(sapply(ontologies, function(x) cons[[x]][,1]), beside=T, names=rep("",3), ylim=c(0,YMAX), space=c(0.2,1),
          border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
          mgp=c(3, 0.75, 0), las=2)

arrows(x0=b,y0=cons.conf.low,y1=cons.conf.high,angle=90,code=3,length=0.05)

## axis labels
mtext("% conserved contacts", side=2, line=2.5,  cex=1)

mtext(c("developmental", "immune", "other"), side=1, at=apply(b, 2, mean), line=0.5, cex=0.8)
mtext("Gene ontology", side=1, line=2, cex=0.8)

## legend & plot label
legend("topright", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],
       fill=dataset.colors[c("Original", "Simulated")], bty='n', 
       inset=c(-0.05, 0), xpd=NA, cex=0.9)

######################################################################################### 

if (ref_sp == "human"){YLIM=c(-0.1,5); tg_sp="mouse"; lab="a"}else{YLIM=c(0, 5); label.enhancers=enhancers; tg_sp="human"; lab="b"}

xlim=c(0.5, length(cons.dist[["all"]][["ENCODE"]]["obs",])+0.5)
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM, xaxs="i", yaxs="i")

for (onto in ontologies){
  points(log(((cons.dist[[onto]][[enh]]["obs",]-cons.dist[[onto]][[enh]]["sim",])/cons.dist[[onto]][[enh]]["sim",])+1),pch=20, col=col.ontologies[onto], cex=0.7)
  
  # Confidence intervals
  segments(x0=1:length(cons.dist.conf.low[[onto]][[enh]]),y0=log(((cons.dist.conf.low[[onto]][[enh]]["obs",]-cons.dist.conf.low[[onto]][[enh]]["sim",])/cons.dist.conf.low[[onto]][[enh]]["sim",])+1),
           x1=1:length(cons.dist.conf.low[[onto]][[enh]]),y1=log(((cons.dist.conf.high[[onto]][[enh]]["obs",]-cons.dist.conf.high[[onto]][[enh]]["sim",])/cons.dist.conf.high[[onto]][[enh]]["sim",])+1),
           col=col.enhancers[enh], lwd=0.3)
}

## axis, legend & plot label
axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0))
mtext("Distance from TSS (Mb)", side=1, line=2)

axis(side=2, mgp=c(3, 0.75, 0), las=2)
mtext("excess of conserved contacts", side=2, line=2.5)

mtext(paste(ref_sp), side=3, line=0.5)

legend("topleft", col=col.ontologies, legend = c("developmental", "immune", "other"), bty='n',pch=20)


dev.off()
