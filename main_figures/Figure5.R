##############################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  
  source("parameters.R")
}

##############################################################################

if(load){
  ref_sp = "human"
  tg=setdiff(c("human", "mouse"), ref_sp)
  
  enhancers=enhancer.datasets[[ref_sp]]
  
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.", ref, ".stats.Rdata", sep=""))
  
  load=FALSE
}

##############################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##############################################################################

#pdf(paste(pathFigures, "Figure5.pdf", sep=""), width=6.85, height=5)

par(mai = c(0.5, 0.5, 0.3, 0.2)) # bottom, left, top, right
mtext.CEX = 0.75

m=matrix(rep(NA, 2*10), nrow=2)
m[1,]=c(rep(1,5), rep(2,5))
m[2,]=c(rep(3,5), rep(4,5))
layout(m)

#################### Fig 5.A - % of conserved contacts #####################

YMAX=40

par(lwd = 1.5)
cons = cons[["all"]]
cons.conf.low = cons.conf.low[["all"]]
cons.conf.high = cons.conf.high[["all"]]

b=barplot(cons, beside=T, names=rep("", dim(cons)[2]), ylim=c(0,YMAX), space=c(0.2,1),
          border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
          mgp=c(3, 0.75, 0), las=2)

arrows(x0=b,y0=cons.conf.low,y1=cons.conf.high,angle=90,code=3,length=0.05)

## axis labels
mtext("% conserved contacts", side=2, line=2.5,  cex=mtext.CEX)

mtext(enh.syn.narrow[enhancers[1:2]],side=1, at=apply(b, 2, mean)[1:2], line=0.5, cex=0.6)
mtext(enh.syn.narrow[enhancers[3:4]],side=1, at=apply(b, 2, mean)[3:4], line=1, cex=0.6)

## legend & plot label
legend("topright", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],
       fill=dataset.colors[c("Original", "Simulated")], bty='n', 
       inset=c(0.05, -0.1), xpd=NA)

mtext("a", side=3, line=1, at=-1.2, font=2, cex=1.05)

############### Fig 5.B - Contact conservation by distance from TSS ##############

YLIM=c(-0.1,5)

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
mtext("distance from promoter region (Mb)", side=1, line=2.5, cex=mtext.CEX)

axis(side=2, mgp=c(3, 0.75, 0), las=2)
mtext("excess of contact conservation", side=2, line=2.5,  cex=mtext.CEX)

abline(h=0, lty=2)
abline(a=0, b=0.1, lty=2) # equation de droite hors log  : y = 0.1*exp(x) ?

legend("topleft", col=col.enhancers, legend = label.enhancers, bty='n',pch=20)
mtext("b", side=3, line=1, at=-3.75, font=2, cex=1.05)

#################  Fig 5.C - Contact conservation by nb cell types #######################

max.nb.cell = 8

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
mtext("number of cell types", side=1, line=2.5, cex=mtext.CEX)

axis(side=2, mgp=c(3, 0.75, 0), las=2)
mtext("% conserved contacts", side=2, line=2.5, cex=mtext.CEX)

mtext("c", side=3, line=1, at=-0.65, font=2, cex=1.05)


#################### Fig 5.D - % of conserved contacts in common cells #####################

YMAX=40
par(lwd = 1.5)

enh = "ENCODE"

b=barplot(cons.common.cell[[enh]], beside=T, names=rep("", dim(cons.common.cell[[enh]])[2]), ylim=c(0,YMAX), space=c(0.2,1),
          border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
          mgp=c(3, 0.75, 0), las=2)

arrows(x0=b,y0=cons.common.cell.conf.low[[enh]],y1=cons.common.cell.conf.high[[enh]],angle=90,code=3,length=0.05)

## axis labels
label.cells = c("ESC", "pre-adipocytes", "B lymphocytes")
mtext("% conserved contacts", side=2, line=2.5, cex=mtext.CEX)
mtext(label.cells, at=apply(b, 2, mean), side=1, line=0.5, cex=mtext.CEX)


mtext("d", side=3, line=1, at=-0.7, font=2, cex=1.05)

####################################################################################

#dev.off()

####################################################################################
