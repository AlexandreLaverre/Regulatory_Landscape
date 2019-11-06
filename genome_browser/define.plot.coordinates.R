######################################################################

define.plot.coordinates<-function(sp.list, focus.genes, gene.coords, annot.baits.TSS, merged.interactions){
  gene.position=list()
  xlim.range=list()

  for(sp in sp.list){
    ## first define gene position
    
    this.gene=focus.genes[sp]
    this.gene.annot=gene.coords[[sp]]

    this.gene.start=this.gene.annot$start[which(this.gene.annot$id==this.gene)]
    this.gene.end=this.gene.annot$end[which(this.gene.annot$id==this.gene)]
    this.chr=this.gene.annot$chr[which(this.gene.annot$id==this.gene)]
    this.strand=this.gene.annot$strand[which(this.gene.annot$id==this.gene)]
    
    this.tss=NA

    if(this.strand=="+"){
      this.tss=this.gene.start
    } else{
      if(this.strand=="-"){
        this.tss=this.gene.end
      } else{
        stop(paste("unknown strand: ",this.strand))
      }
    }

    gene.position[[sp]]=list("chr"=this.chr, "strand"=this.strand, "start"=this.gene.start, "end"=this.gene.end, "tss"=this.tss)

    ## bait position

    this.bait.annot=annot.baits.TSS[[sp]]
    this.bait.id=this.bait.annot[which(this.bait.annot$gene_ID==this.gene), "bait_ID"]

    ## merged interactions

    this.merged.int=merged.interactions[[sp]]
    this.merged.int=this.merged.int[which(this.merged.int[,"bait_ID"]%in%this.bait.id),]

    this.int.range=range(c(this.merged.int$bait_start, this.merged.int$bait_end, this.merged.int$start, this.merged.int$end))
    this.total.range=range(c(this.int.range, this.gene.start, this.gene.end))

    if(this.strand=="+"){
      size.5prime=this.tss-this.total.range[1]
      size.3prime=this.total.range[2]-this.tss
    } else{
      size.5prime=this.total.range[2]-this.tss
      size.3prime=this.tss-this.total.range[1]
    }
     
    xlim.range[[sp]]=c("size.5prime"=size.5prime, "size.3prime"=size.3prime)
  }

  max.size5=max(unlist(lapply(xlim.range, function(x) x[["size.5prime"]])))
  max.size3=max(unlist(lapply(xlim.range, function(x) x[["size.3prime"]])))

  final.coordinates=list()

  for(sp in sp.list){
    this.chr=gene.position[[sp]][["chr"]]
    this.strand=gene.position[[sp]][["strand"]]
    this.tss=gene.position[[sp]][["tss"]]

    if(this.strand=="+"){
      this.start=this.tss-max.size5
      this.end=this.tss+max.size3
    } else{
      ## start > end
      this.start=this.tss+max.size5
      this.end=this.tss-max.size3 
    }

    final.coordinates[[sp]]=list("chr"=this.chr, "strand"=this.strand, "tss"=this.tss, "start"=this.start, "end"=this.end)
  }

  return(final.coordinates)

}

######################################################################
