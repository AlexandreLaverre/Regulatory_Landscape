##############################################################################
setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")


source("parameters.R")

nb = 1
mtext.CEX = 0.9

pdf(paste(pathFigures, "SupplementaryFigure28.pdf", sep=""), width=7, height=5)

par(mai = c(0.5, 0.45, 0.5, 0.1)) # bottom, left, top, right
par(lwd = 1.5)

par(mfrow=c(2,3))

############### Supplementary Fig 28 - Contact conservation by distance from TSS ##############

for (ref_sp in c("human", "mouse")){
  
  enhancers=enhancer.datasets[[ref_sp]]
  load(paste(pathFigures, "RData/data.contact.conservation.", ref_sp, ".Rdata", sep=""))
  
  if (ref_sp == "human"){YMAX=40; tg_sp="mouse"; enhancers=enhancers[-1]}else{YMAX=25; tg_sp="human"}
  
  for (enh in enhancers){
    b=barplot(cons.common.cell[[enh]], beside=T, names=rep("", dim(cons.common.cell[[enh]])[2]), ylim=c(0,YMAX), space=c(0.2,1),
              border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
              mgp=c(3, 0.75, 0), las=2, xaxt='n')
    
    arrows(x0=b,y0=cons.common.cell.conf.low[[enh]],y1=cons.common.cell.conf.high[[enh]],angle=90,code=3,length=0.05)
    
    ## axis labels
    label.cells = c("ESC", "Pre-\nadipocytes", "Bcell")
    mtext(label.cells, at=apply(b, 2, mean), side=1, line=1.5, cex=0.8)
    
    if (nb == 1 | nb == 4){axis(side=2, mgp=c(3, 0.75, 0), las=2)
      mtext("% of conserved contacts", side=2, line=2, cex=mtext.CEX)}
    
    mtext(paste(enh, "\n", ref_sp, "vs.", tg_sp), side=3, line=0.5, cex=0.8)
    
    mtext(letters[nb], side=3, line=1, at=0.1, font=2, cex=mtext.CEX)
    nb = nb +1
  }
  
}

plot.new()
legend("center", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],
       fill=dataset.colors[c("Original", "Simulated")], bty='n', inset=c(0.05, -0.1), xpd=NA, cex=1.5)
############################################# 
dev.off()
