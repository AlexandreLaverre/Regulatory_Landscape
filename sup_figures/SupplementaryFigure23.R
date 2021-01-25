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
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.", ref, ".stats.Rdata", sep=""))
  
  enh="ENCODE"
  
  load=FALSE
}

##############################################################################

if(prepare){
  
  sampleinfo.ref=sampleinfo[[ref]]
  sampleinfo.tg=sampleinfo[[tg]]
  
  rownames(sampleinfo.ref)=sampleinfo.ref$Sample.ID
  rownames(sampleinfo.tg)=sampleinfo.tg$Sample.ID
  
  ## excess of contact conservation by sample
  
  mat.cons.obs=100*matrix.cons.obs[[enh]]
  mat.cons.sim=100*matrix.cons.sim[[enh]]
  
  diffmat=mat.cons.obs-mat.cons.sim
  
  s.ref=rep(rownames(diffmat), dim(diffmat)[2])
  s.tg=rep(colnames(diffmat), each=dim(diffmat)[1])
  
  diffmat.table=data.frame("diff"=as.numeric(diffmat), "sample.ref"=s.ref, "sample.tg"=s.tg, "cell.ref"=sampleinfo.ref[s.ref,"Broad.cell.type.or.tissue"], "cell.tg"=sampleinfo.tg[s.tg,"Broad.cell.type.or.tissue"])
  
  diffmat=t(diffmat)
  
  ## conservation by distance
  cons.dist = cons.dist[["all"]]
  cons.dist.conf.low = cons.dist.conf.low[["all"]]
  cons.dist.conf.high = cons.dist.conf.high[["all"]]
  
  
  prepare=FALSE
}

##############################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##############################################################################

pdf(paste(pathFigures, "SupplementaryFigure23.pdf", sep=""), width=6.85, height=6.5)

m=matrix(rep(NA, 9*20), nrow=9)

for(i in 1:5){
  m[i,]=rep(1,20)
}

for(i in 6:9){
  m[i,]=c(rep(2,10), rep(3, 10))
}



layout(m)

mtext.CEX=0.7

##############################################################################

## heatmap for excess of contact conservation by cell type

zmax=max(abs(as.numeric(diffmat)))
zlim=c(-zmax, zmax)   ## range(as.numeric(diffmat))

par(mar=c(9, 12.1, 2, 1.1))

nbcol=50

col.heatmap=colorRampPalette(c("darkred", "white", "navy"), space="rgb")(nbcol)

image(diffmat, zlim=zlim, axes=F, col=col.heatmap)

xpos=seq(from=0, to=1, length=dim(diffmat)[1])
names(xpos)=rownames(diffmat)

ypos=seq(from=0, to=1, length=dim(diffmat)[2])
names(ypos)=colnames(diffmat)

smally=diff(ypos)[1]/2
smallx=diff(xpos)[1]/2

tinyx=diff(xpos)[1]/10
tinyy=diff(ypos)[1]/10

for(c in unique(sampleinfo.ref$Broad.cell.type.or.tissue)){
  s=sampleinfo.ref$Sample.ID[which(sampleinfo.ref$Broad.cell.type.or.tissue==c)]
  this.y=ypos[s]
  mtext(syn.celltypes[c], side=2, line=0.75, cex=0.7, las=2, at=mean(this.y))
  
  segments(-0.03, min(this.y)-tinyy,  -0.03, max(this.y)+tinyy, xpd=NA)
}

par(xpd = TRUE)
for(c in unique(sampleinfo.tg$Broad.cell.type.or.tissue)){
  s=sampleinfo.tg$Sample.ID[which(sampleinfo.tg$Broad.cell.type.or.tissue==c)]
  this.x=xpos[s]
  text(mean(this.x), y =-0.1, syn.celltypes[c], srt = 45, cex=1, adj=1)
  
  segments(min(this.x)-tinyx, -0.07,  max(this.x)+tinyx,  -0.07, xpd=NA)
}

par(xpd = FALSE)

for(c in c("pre-adipocytes", "embryonic stem cells", "B lymphocytes")){
  sr=sampleinfo.ref$Sample.ID[which(sampleinfo.ref$Broad.cell.type.or.tissue==c)]
  this.y=ypos[sr]
  
  st=sampleinfo.tg$Sample.ID[which(sampleinfo.tg$Broad.cell.type.or.tissue==c)]
  this.x=xpos[st]
  
  rect(min(this.x)-smallx, min(this.y)-smally, max(this.x)+smallx, max(this.y)+smally, border="gold", xpd=NA)
}

## legend for the heatmap

startx=seq(from=-0.32, to=-0.12, length=nbcol)
step=diff(startx)[1]
endx=startx+step

starty=rep(-0.15, 1)
endy=rep(-0.1, 1)
tty=(endy-starty)/4.5
sy=(endy-starty)/1.25

rect(startx, starty, endx, endy, col=col.heatmap, border=col.heatmap, xpd=NA)

segments(min(startx), starty-tty, max(endx), starty-tty, xpd=NA)
segments(min(startx), starty-sy,  min(startx), starty-tty, xpd=NA)
segments(max(endx), starty-sy,  max(endx), starty-tty, xpd=NA)
segments(mean(c(startx,endx)), starty-sy,  mean(c(startx,endx)), starty-tty, xpd=NA)

text(round(zlim[1], digits=0), x=min(startx), y=starty-1.75*sy, xpd=NA, cex=0.95)
text(round(zlim[2], digits=0), x=max(endx), y=starty-1.75*sy, xpd=NA, cex=0.95)
text("0.0", x=mean(c(startx,endx)), y=starty-1.75*sy, xpd=NA, cex=0.95)

text("% contact conservation", x=mean(c(startx,endx)), y=starty-3.2*sy, xpd=NA, cex=1.07)
text("(observed-simulated)", x=mean(c(startx,endx)), y=starty-4.7*sy, xpd=NA, cex=1.07)

mtext("a", side=3, font=2, line=0.75, at=-0.3, cex=1)

##############################################################################

## conservation as a number of cell types

max.nb.cell = 6

ylim=c(10, 60)
xlim=c(0.5, max.nb.cell+0.5)
xpos=seq(1, max.nb.cell, 1)
names(xpos) = 1:max.nb.cell

par(mar=c(3.5, 4.5, 2, 2.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(nb_cell in 1:max.nb.cell){
  x=xpos[nb_cell]
  
  points(x, cons.nb.cell[[enh]]["obs",nb_cell], pch=20, col=dataset.colors["Original"])
  segments(x, cons.nb.cell.conf.low[[enh]]["obs", nb_cell], x, cons.nb.cell.conf.high[[enh]]["obs", nb_cell], col=dataset.colors["Original"])
  
  points(x, cons.nb.cell[[enh]]["sim",nb_cell], pch=20, col=dataset.colors["Simulated"])
  segments(x, cons.nb.cell.conf.low[[enh]]["sim", nb_cell], x, cons.nb.cell.conf.high[[enh]]["sim", nb_cell], col=dataset.colors["Simulated"])
}

abline(v=xpos[1:max.nb.cell-1]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", max.nb.cell), mgp=c(3, 0.5, 0))
mtext(colnames(cons.nb.cell[[enh]]), at=xpos, side=1, line=0.5, cex=mtext.CEX)
mtext("number of cell types", side=1, line=2, cex=mtext.CEX)

axis(side=2, mgp=c(3, 0.75, 0), las=2)
mtext("% conserved contacts", side=2, line=2.5, cex=mtext.CEX)

legend("topleft", col=dataset.colors, legend = c("PCHi-C data", "simulated data"), box.col="white", bg="white", pch=20, inset=c(0.01, 0))

mtext("b", side=3, line=1.75, at=-0.3, font=2, cex=1)

##############################################################################

YLIM=c(-1,35)

class_leg <- c("0",  "0.5",  "1", "1.5", "2")

par(lwd = 0.7)

xlim=c(0.5, length(cons.dist[["ENCODE"]]["obs",])+0.5)

par(mar=c(3.5, 4.5, 2, 2.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM,  yaxs="i")
nbclass=length(cons.dist.conf.low[[enh]]["obs",])

points(cons.dist[[enh]]["obs",],pch=20, col=dataset.colors["Original"])
segments(1:nbclass, cons.dist.conf.low[[enh]]["obs",], 1:nbclass, cons.dist.conf.high[[enh]]["obs",], col=dataset.colors["Original"])

points(cons.dist[[enh]]["sim",],pch=20, col=dataset.colors["Simulated"])
segments(1:nbclass, cons.dist.conf.low[[enh]]["sim",], 1:nbclass, cons.dist.conf.high[[enh]]["sim",], col=dataset.colors["Simulated"])

## axis, legend & plot label
axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0))
mtext("distance from promoter region (Mb)", side=1, line=2, cex=mtext.CEX)

axis(side=2, mgp=c(3, 0.75, 0), las=2)
mtext("% conserved contacts", side=2, line=2.5,  cex=mtext.CEX)

mtext("c", side=3, line=1.75, at=-7.5, font=2, cex=1.05)

##############################################################################

dev.off()

##############################################################################
##############################################################################
