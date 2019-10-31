######################################################################

objects=ls()

if(!("load"%in%objects)){
  load=TRUE
  process=TRUE
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

source("plot.annotations.R")

######################################################################

if(process==TRUE){
  
  sp="human"
  samples= names(interactions[[sp]])
  nbsamples=length(samples)
  gene="ENSG00000164690"
  
  gene.annot=gene.coords[[sp]]
  this.gene.start=gene.annot$start[which(gene.annot$id==gene)]
  this.gene.end=gene.annot$end[which(gene.annot$id==gene)]
  this.chr=gene.annot$chr[which(gene.annot$id==gene)]

  exon.annot=exon.coords[[sp]]
  
  bait.annot=annot.baits.TSS[[sp]]
  bait.thisgene=bait.annot[which(bait.annot$gene_ID==gene),]

  ## interactions by sample, before merging fragments
  
  interactions.thisgene=data.frame()
  
  for(sample in samples){
    this.int=interactions[[sp]][[sample]]
    this.int$bait_ID=paste(this.int$bait_chr, this.int$bait_start, this.int$bait_end, sep=",")
    this.int=this.int[which(this.int$bait_ID%in%bait.thisgene$bait_ID),]
    
    if(dim(this.int)[1]>0){
      this.int$sample=rep(sample, dim(this.int)[1])
      
      if(dim(interactions.thisgene)[1]>0){
        interactions.thisgene=rbind(interactions.thisgene, this.int, stringsAsFactors=F)
      } else{
        interactions.thisgene=this.int
      }
    }
  }

  ## merged interactions

  this.merged.int=merged.interactions[[sp]][which(merged.interactions[[sp]][,"bait_ID"]%in%bait.thisgene$bait_ID),]
  
  
  xstart=min(c(this.gene.start, interactions.thisgene$otherEnd_start))
  xend=max(c(this.gene.end, interactions.thisgene$otherEnd_end))


  ## enhancer coordinates

  this.enhancers=enhancer.coords[[sp]]
  
  process=TRUE
}

######################################################################

panelheight.enhancers=2
panelheight.annotations=3
panelheight.interactions=4

nbpanels=panelheight.annotations+panelheight.enhancers+panelheight.interactions

figheight=nbpanels*0.25

pdf(file=paste("figures/", gene, "_",xstart,"_", xend, ".pdf",sep=""), width=8, height=figheight)

######################################################################

m=matrix(rep(NA, nbpanels), nrow=nbpanels)

for(i in 1:panelheight.annotations){
  m[i,]=1
}

for(i in (panelheight.annotations+1):(panelheight.annotations+panelheight.enhancers)){
  m[i,]=2
}


for(i in (panelheight.annotations+panelheight.enhancers+1):nbpanels){
  m[i,]=3
}

layout(m)

######################################################################

## plot annotations

par(mar=c(1.5, 2.1, 0.25, 1.1))
plot.annotations.genes(gene.annot, gene.biotypes=c("protein_coding", "lincRNA", "lncRNA", "processed_transcript"), focus.gene=gene, xlim=c(xstart, xend), axis=T, axisunit="Mb")

######################################################################

## plot enhancers

plot.enhancers(this.enhancers, chr=this.chr, xlim=c(xstart, xend), col="gray40")

######################################################################

## plot interactions



######################################################################

## for(sample in samples){
##   plot(1, type="n", xlim=c(xstart, xend), ylim=c(0, 1), xaxs="i", yaxs="i", xlab="", axes=F)
  
##   this.int=interactions.thisgene[which(interactions.thisgene$sample==sample),]

##   if(dim(this.int)[1]>0){
##     rect(this.int$otherEnd_start, 0.15, this.int$otherEnd_end, 0.95, col="red", border="red")
##   }
  
## }

######################################################################

dev.off()

######################################################################
