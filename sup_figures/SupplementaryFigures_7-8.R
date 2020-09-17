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

  load=FALSE
}

##########################################################################

if(prepare){

  fignb=c(7, 8)
  names(fignb)=c("human", "mouse")

  prepare=FALSE
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##########################################################################

for(sp in c("human", "mouse")){
  pdf(paste(pathFigures, "SupplementaryFigure",fignb[sp], ".pdf", sep=""), width=4.49, height=3.5)

  ## layout
  
  m=matrix(rep(NA,11*23), nrow=11)
  for(i in 1:10){
    m[i,]=c(rep(1,5), rep(2, 15), rep(3, 3))
  }
  for(i in 11){
    m[i,]=c(rep(4, 23))
  }
  layout(m)

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
  par(mar=c(0.8,1,1.5,0.1))
  plot(tree, direction="rightwards", show.tip.label=FALSE)

  ## matrix obs-sim

  par(mar=c(1,0,2,0))
  image(mat.diff, axes=F, col=terrain.colors(50), zlim=c(0, 100))

  ## colors by cell type

  par(mar=c(1, 0, 2, 0))
  plot(1, type="n", xlab="", ylab="", main="", axes=F)

  this.samples=rownames(mat.diff)
  this.cells=info[this.samples, 
  
  dev.off()
}

##########################################################################
