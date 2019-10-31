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

  exon.annot=exon.coords[[sp]]
  
  bait.annot=annot.baits.TSS[[sp]]
  bait.thisgene=bait.annot[which(bait.annot$gene_ID==gene),]
  
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
  
  xstart=min(c(this.gene.start, interactions.thisgene$otherEnd_start))
  xend=max(c(this.gene.end, interactions.thisgene$otherEnd_end))

  process=TRUE
}

######################################################################

figheight=(nbsamples+1)*0.25

pdf(file=paste("figures/", gene, "_",xstart,"_", xend, ".pdf",sep=""), width=8, height=figheight)

######################################################################

m=matrix(rep(NA, nbsamples+2), nrow=nbsamples+2)

for(i in 1:2){
  m[i,]=1
}

for(i in 3:(nbsamples+2)){
  m[i,]=i-1
}

layout(m)

######################################################################

par(mar=c(1.5, 2.1, 0.25, 1.1))
plot.annotations(gene.annot, focus.gene=gene, xlim=c(xstart, xend))

######################################################################

for(sample in samples){
  plot(1, type="n", xlim=c(xstart, xend), ylim=c(0, 1), xaxs="i", yaxs="i", xlab="", axes=F)
  
  this.int=interactions.thisgene[which(interactions.thisgene$sample==sample),]

  if(dim(this.int)[1]>0){
    rect(this.int$otherEnd_start, 0.15, this.int$otherEnd_end, 0.95, col="red", border="red")
  }
  
}

######################################################################

dev.off()

######################################################################
