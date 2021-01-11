##############################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  
  source("../main_figures/parameters.R")
}

##############################################################################

if(load){
  ref_sp = "mouse"
  tg=setdiff(c("human", "mouse"), ref_sp)
  
  enhancers=enhancer.datasets[[ref_sp]]
  
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.", ref_sp, ".stats.Rdata", sep=""))
  
  load=FALSE
}

##############################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##############################################################################

#pdf(paste(pathFigures, "SupplementaryFigure24.pdf", sep=""), width=6.85, height=5)

par(mai = c(0.5, 0.5, 0.3, 0.2)) # bottom, left, top, right
mtext.CEX = 0.75

m=matrix(rep(NA, 2*10), nrow=2)
m[1,]=c(rep(1,5), rep(2,5))
m[2,]=c(rep(3,5), rep(4,5))
layout(m)

#################### Fig 5.A - % of conserved contacts #####################

YMAX=45

par(lwd = 1.5)

cons = cons[["all"]]
cons.conf.low = cons.conf.low[["all"]]
cons.conf.high = cons.conf.high[["all"]]


b=barplot(cons, beside=T, names=rep("", dim(cons)[2]), ylim=c(0,YMAX), space=c(0.2,1),
          border=dataset.colors[c("Original", "Simulated")],  col=c("yellow", "cadetblue1", dataset.colors[c("Original", "Simulated")]),
          mgp=c(3, 0.75, 0), las=2)

arrows(x0=b,y0=cons.conf.low,y1=cons.conf.high,angle=90,code=3,length=0.05)

## axis labels
mtext("% of conserved contacts", side=2, line=2.5,  cex=mtext.CEX)
mtext(c("ENCODE", "FANTOM5"), side=1, at=apply(b[,1:2], 2, mean), line=0.5, cex=0.6)
if (ref_sp == "human"){mtext(c("FOCS\nGRO-seq", "RoadMap\nEpigenomics"), side=1, at=apply(b[,3:4], 2, mean), line=1.2, cex=0.6)}

## legend & plot label
legend("topright", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],
       fill=dataset.colors[c("Original", "Simulated")], bty='n', 
       inset=c(0.05, -0.1), xpd=NA)

mtext("a", side=3, line=1, at=0.1, font=2, cex=1.05)

############### Fig 5.B - Contact conservation by distance from TSS ##############
YLIM=c(-0.2, 8)
label.enhancers=enhancers
class_leg <- c("0",  "0.5",  "1", "1.5", "2")

par(lwd = 0.7)
cons.dist = cons.dist[["all"]]
cons.dist.conf.low = cons.dist.conf.low[["all"]]
cons.dist.conf.high = cons.dist.conf.high[["all"]]

xlim=c(0.5, length(cons.dist[["ENCODE"]]["obs",])+0.5)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM, xaxs="i", yaxs="i")

for (enh in enhancer.datasets[[ref_sp]]){
  
  points(log(((cons.dist[[enh]]["obs",]-cons.dist[[enh]]["sim",])/cons.dist[[enh]]["sim",])+1),pch=20, col=col.enhancers[enh])
  
  # Confidence intervals
  segments(x0=1:length(cons.dist.conf.low[[enh]]),y0=log(((cons.dist.conf.low[[enh]]["obs",]-cons.dist.conf.low[[enh]]["sim",])/cons.dist.conf.low[[enh]]["sim",])+1),
           x1=1:length(cons.dist.conf.low[[enh]]),y1=log(((cons.dist.conf.high[[enh]]["obs",]-cons.dist.conf.high[[enh]]["sim",])/cons.dist.conf.high[[enh]]["sim",])+1),
           col=col.enhancers[enh], lwd=0.3)
}

## axis, legend & plot label
axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0))
mtext("Distance from TSS (Mb)", side=1, line=2.5, cex=mtext.CEX)

axis(side=2, mgp=c(3, 0.75, 0), las=2)
mtext("Excess of contact conservation", side=2, line=2.5,  cex=mtext.CEX)
abline(h=0, lty=2)

legend("topleft", col=col.enhancers, legend = label.enhancers, bty='n',pch=20)
mtext("b", side=3, line=1, at=0.1, font=2, cex=1.05)

#################  Fig 5.C - Contact conservation by nb cell types #######################
if (ref_sp == "human"){max.nb.cell = 8}else{max.nb.cell = 6}

ylim=c(10, 70)
xlim=c(0.5, max.nb.cell+0.5)
xpos=seq(1, max.nb.cell, 1)
names(xpos) = 1:max.nb.cell

smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[ref_sp]]

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[ref_sp]]){
  for(nb_cell in 1:max.nb.cell){
    
    x=xpos[nb_cell]+smallx[enh]
    points(x, cons.nb.cell[[enh]]["obs",nb_cell], pch=20, col=col.enhancers[enh])
    segments(x, cons.nb.cell.conf.low[[enh]]["obs", nb_cell], x, cons.nb.cell.conf.high[[enh]]["obs", nb_cell], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:max.nb.cell-1]+0.5, lty=3, col="gray40")
axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", max.nb.cell))
mtext(colnames(cons.nb.cell[[enh]]), at=xpos, side=1, line=1, cex=mtext.CEX)
mtext("Number of cell types", side=1, line=2.5, cex=mtext.CEX)

axis(side=2, mgp=c(3, 0.75, 0), las=2)
mtext("Conserved contact (%)", side=2, line=2.5, cex=mtext.CEX)

legend("topleft", legend=label.enhancers, pch=20,
       col=col.enhancers, bty="n", cex=1)

mtext("c", side=3, line=1, at=0.1, font=2, cex=1.05)


#################### Fig 5.D - % of conserved contacts in common cells #####################
if (ref_sp == "human"){YMAX=40}else{YMAX=25}
par(lwd = 1.5)

enh = "ENCODE"

b=barplot(cons.common.cell[[enh]], beside=T, names=rep("", dim(cons.common.cell[[enh]])[2]), ylim=c(0,YMAX), space=c(0.2,1),
          border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
          mgp=c(3, 0.75, 0), las=2)

arrows(x0=b,y0=cons.common.cell.conf.low[[enh]],y1=cons.common.cell.conf.high[[enh]],angle=90,code=3,length=0.05)

## axis labels
label.cells = c("ESC", "Preadipocytes", "Bcell")
mtext("% of conserved contacts", side=2, line=2.5, cex=mtext.CEX)
mtext(label.cells, at=apply(b, 2, mean), side=1, line=0.5, cex=mtext.CEX)

## legend & plot label
legend("topright", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],
       fill=dataset.colors[c("Original", "Simulated")], bty='n', 
       inset=c(0.05, -0.1), xpd=NA)

mtext("d", side=3, line=1, at=0.1, font=2, cex=1.05)

####################################################################################

dev.off()

####################################################################################
