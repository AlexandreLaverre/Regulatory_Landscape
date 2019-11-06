######################################################################

objects=ls()

if(!("load"%in%objects)){
  load=TRUE
  process=TRUE
  
  source("plot.annotations.R")
  source("plot.enhancers.R")
  source("plot.interactions.R")
}

######################################################################

if(load==TRUE){
  load("RData/data.annotations.RData")
  load("RData/data.interactions.per.sample.RData")
  load("RData/data.interactions.annotations.RData")
  load("RData/data.merged.interactions.RData")
  load("RData/data.enhancers.RData")

  load=FALSE
}
######################################################################

if(process==TRUE){

  genes=c("ENSG00000164690", "ENSMUSG00000002633")
  names(genes)=c("human", "mouse")

  plot.coords=define.plot.coordinates(sp.list=c("human", "mouse"), focus.genes=genes, gene.coords=gene.coords, annot.baits.TSS=annot.baits.TSS, merged.interactions=merged.interactions)


  print(plot.coords)
  
  ## sp="human"
  ## samples=names(interactions[[sp]])
  ## nbsamples=length(samples)
  ## gene=
  
  ## gene.annot=gene.coords[[sp]]
  ## this.gene.start=gene.annot$start[which(gene.annot$id==gene)]
  ## this.gene.end=gene.annot$end[which(gene.annot$id==gene)]
  ## this.chr=gene.annot$chr[which(gene.annot$id==gene)]

  ## exon.annot=exon.coords[[sp]]
   
  ## bait.annot=annot.baits.TSS[[sp]]
  ## bait.thisgene=bait.annot[which(bait.annot$gene_ID==gene),]

  ## ## interactions by sample, before merging fragments
  
  ## interactions.thisgene=data.frame()
  
  ## for(sample in samples){
  ##   this.int=interactions[[sp]][[sample]]
  ##   this.int$bait_ID=paste(this.int$bait_chr, this.int$bait_start, this.int$bait_end, sep=",")
  ##   this.int=this.int[which(this.int$bait_ID%in%bait.thisgene$bait_ID),]
    
  ##   if(dim(this.int)[1]>0){
  ##     this.int$sample=rep(sample, dim(this.int)[1])
      
  ##     if(dim(interactions.thisgene)[1]>0){
  ##       interactions.thisgene=rbind(interactions.thisgene, this.int, stringsAsFactors=F)
  ##     } else{
  ##       interactions.thisgene=this.int
  ##     }
  ##   }
  ## }

  ## ## merged interactions

  ## this.merged.int=merged.interactions[[sp]][which(merged.interactions[[sp]][,"bait_ID"]%in%bait.thisgene$bait_ID),]
  
  
  ## xstart=min(c(this.gene.start, interactions.thisgene$otherEnd_start))
  ## xend=max(c(this.gene.end, interactions.thisgene$otherEnd_end))


  ## ## enhancer coordinates

  ## this.enhancers=enhancer.coords[[sp]]
  
  process=FALSE
}

######################################################################

## panelheight.enhancers=4
## panelheight.annotations=6
## panelheight.interactions=2

## nbpanels=panelheight.annotations+panelheight.enhancers+panelheight.interactions

## figheight=nbpanels*0.15

## pdf(file=paste("figures/", gene, "_",xstart,"_", xend, ".pdf",sep=""), width=8, height=figheight)

## ######################################################################

## m=matrix(rep(NA, nbpanels), nrow=nbpanels)

## for(i in 1:panelheight.annotations){
##   m[i,]=1
## }

## for(i in (panelheight.annotations+1):(panelheight.annotations+panelheight.enhancers)){
##   m[i,]=2
## }


## for(i in (panelheight.annotations+panelheight.enhancers+1):nbpanels){
##   m[i,]=3
## }

## layout(m)

## ######################################################################

## ## plot annotations

## par(mar=c(1.5, 2.1, 0.25, 1.1))

## plot.annotations.genes(gene.annot, gene.biotypes=c("protein_coding"), focus.gene=gene, xlim=c(xstart, xend), axis=T, axisunit="Mb")

## ######################################################################

## ## plot enhancers

## par(mar=c(0.5, 2.1, 0.5, 1.1))
## plot.enhancers(this.enhancers, chr=this.chr, xlim=c(xstart, xend), col="gray40", separate=TRUE)

## ######################################################################

## ## plot interactions

## par(mar=c(0.5, 2.1, 0.5, 1.1))
## plot.interactions(this.merged.int, xlim=c(xstart, xend), focus.bait=bait.thisgene$baitID, col.contact="navy", col.otherbait="red")
              
## ######################################################################

## ## for(sample in samples){
## ##   plot(1, type="n", xlim=c(xstart, xend), ylim=c(0, 1), xaxs="i", yaxs="i", xlab="", axes=F)
  
## ##   this.int=interactions.thisgene[which(interactions.thisgene$sample==sample),]

## ##   if(dim(this.int)[1]>0){
## ##     rect(this.int$otherEnd_start, 0.15, this.int$otherEnd_end, 0.95, col="red", border="red")
## ##   }
  
## ## }

## ######################################################################

## dev.off()

## ######################################################################
