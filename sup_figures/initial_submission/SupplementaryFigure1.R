#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

##########################################################################

if(load){
  library(ape)
  
  load(paste(pathFigures, "RData/data.sample.clustering.RData", sep=""))
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.AFC.RData", sep=""))

  load=FALSE
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "SupplementaryFigure1.pdf", sep=""), width=6.85, height=5.5)

## layout
m=matrix(rep(NA,20*40), nrow=20)
for(i in 1:9){
  m[i,]=c(rep(1,3), rep(2, 20), rep(3, 17))
}
m[10,]=c(rep(1,3), rep(2, 20), rep(4, 17))

for(i in 11:19){
  m[i,]=c(rep(5,3), rep(6, 20), rep(7, 17))
}
m[20,]=c(rep(5,3), rep(6, 20), rep(8, 17))

layout(m)

##########################################################################
fig = 1

for(sp in c("human", "mouse")){
  ## get data for this species
  info=sampleinfo[[sp]]
  rownames(info)=info$Sample.ID
  
  tree=as.phylo(sample.clustering[[sp]][["hclust.alldist"]])
  sample.order=sample.clustering[[sp]][["sample.order.alldist"]]
  mat.obs=sample.clustering[[sp]][["mat.alldist.obs"]]
  mat.sim=sample.clustering[[sp]][["mat.alldist.sim"]]
  mat.diff=as.matrix(mat.obs-mat.sim)
  diag(mat.diff)=rep(NA, length(sample.order))
  
  mat.diff=mat.diff[sample.order, sample.order]

  mat.diff=100*mat.diff ## percentage instead of fraction
  
  ## tree
  par(mar=c(0.75,1,1.35,0.1))
  plot(tree, direction="rightwards", show.tip.label=FALSE)

  ## plot label
  mtext(letters[fig], side=3, at=0, font=2, line=0, cex=0.95)
  fig = fig+1
  
  ## matrix obs-sim

  par(mar=c(1,0,1.75,8.25))
  image(mat.diff, axes=F, col=terrain.colors(50), zlim=c(0, 100))

  ## colors by cell type

  this.samples=rownames(mat.diff)
  this.cells=info[this.samples, "Broad.cell.type.or.tissue"]
  names(this.cells)=this.samples

  ypos=seq(from=0, to=1, length=dim(mat.diff)[1])
  ywidth=diff(ypos)[1]

  for(c in unique(this.cells)){
    all.ypos=ypos[which(this.cells==c)]

    if(diff(range(which(this.cells==c)))==(length(all.ypos)-1)){
      ## perfect clustering
      segments(1+ywidth*0.75, min(all.ypos)-ywidth/3, 1+ywidth*0.75, max(all.ypos)+ywidth/3, xpd=NA)
      mtext(syn.celltypes[c], side=4, line=0.75, las=2, cex=0.6, at=mean(all.ypos), col=col.celltypes[c])
    } else{
      segments(1+ywidth*0.75, all.ypos-ywidth/3, 1+ywidth*0.75, all.ypos+ywidth/3, xpd=NA)
      mtext(syn.celltypes[c], side=4, line=0.75, las=2, cex=0.6, at=all.ypos, col=col.celltypes[c])
    }
  }

  mtext(paste("hierarchical clustering", sp, sep=", "), side=3, line=0.5, cex=0.7)

  ## AFC plot

  afc=data.AFC[[sp]][["AFC"]]
  explained=round(100*afc$eig/sum(afc$eig), digits=1)

  par(mar=c(4.1, 5.5, 2, 1.5))
  plot(afc$li[,1], afc$li[,2], pch=20, col=col.celltypes[this.cells[rownames(afc$li)]], xlab="", ylab="", axes=F, cex=1.25)
  box()

  xlim=range(afc$li[,1])

  axis(side=1, mgp=c(3, 0.5, 0), cex.axis=0.85)
  axis(side=2, mgp=c(3, 0.5, 0), cex.axis=0.85)

  mtext(paste("correspondence analysis", sp, sep=", "), side=3, line=0.5, cex=0.7)
  mtext(paste("axis 1 (", explained[1],"% explained variance)",sep=""), side=1, line=1.75, cex=0.7)
  mtext(paste("axis 2 (", explained[2],"% explained variance)",sep=""), side=2, line=1.75, cex=0.7)

  ## plot labels
  mtext(letters[fig], side=3, at=xlim[1]-diff(xlim)/5.3, font=2, line=0.5, cex=0.95)
  fig = fig+1
  
  ## legend for the heatmap

  if(sp=="mouse"){
    par(mar=c(1.45,1.1,0.0,13.1))
    z=seq(0, 100, length = 50)
    zlim=c(0,100)
    xax=c(0, 25, 50, 75, 100)
    xax=xax[which(xax>=min(z) & xax<=max(z))]
    image(x=z, z = matrix(z, ncol = 1), col = terrain.colors(50), zlim=zlim, xlim=range(xax)+c(-2,2), xaxt="n" ,yaxt="n")
    
    par(tck=-0.75)
    axis(side=1, at = xax, labels = xax, cex.axis=0.85, mgp=c(3,0.4,0))
    
    mtext("% shared interactions", side=4, las=2, at=1.15, cex=0.65, line=1)
    mtext("(observed-simulated)", side=4, las=2, at=-1.55, cex=0.65, line=1)
    par(tck=NA)
  } else{
    ## empty plot
    par(mar=c(1.35,1.1,0.0,13.1))
    plot(1, type="n", xlab="", ylab="", main="", axes=F)
  }

}

##########################################################################

dev.off()
