#####################################################################

plot.annotations.genes<-function(gene.coords, focus.gene, gene.biotypes="all", xlim, col.focus="navy", col.other="gray60", axis=T){

  if(!focus.gene%in%gene.coords$id){
    stop(paste("cannot find", focus.gene, "in annotations"))
  }
  ## extract chr
  this.chr=gene.coords$chr[which(gene.coords$id==focus.gene)]
  this.genes=gene.coords[which(gene.coords$chr==this.chr & ((gene.coords$start>=xlim[1] & gene.coords$start<=xlim[2]) | (gene.coords$end>=xlim[1] & gene.coords$end<=xlim[2]) | (gene.coords$start<=xlim[1] & gene.coords$end>=xlim[2]))),]

  if(gene.biotypes!="all"){
    this.genes=this.genes[which(this.genes$biotype%in%gene.biotypes),]
  }
  
  if(!focus.gene%in%this.genes$id){
   stop(paste(focus.gene, "does not overlap with region")) 
  }

  this.genes$color=col.other
  this.genes$color[which(this.genes$id==focus.gene)]=col.focus
  
  print(paste(dim(this.genes)[1], "genes"))
  
  ylim=c(-2,2)
  ypos.bystrand=c(1,-1)
  names(ypos.bystrand)=c("1", "-1")
  height=0.2

  arrowheight=0.25
  arrowsize=diff(xlim)/50
  
  plot(1, type="n", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, xaxs="i", yaxs="i") ## empty plot
  
  ## abline(h=0, col="gray40")

  for(g in this.genes$id){
    this.strand=as.character(this.genes$strand[which(this.genes$id==g)])
    this.col=this.genes$color[which(this.genes$id==g)]
    this.start=this.genes$start[which(this.genes$id==g)]
    this.end=this.genes$end[which(this.genes$id==g)]
    this.name=this.genes$name[which(this.genes$id==g)]

    ypos=ypos.bystrand[this.strand]

    if(this.start<xlim[1]){
      this.start=xlim[1]
    }

    if(this.end>xlim[2]){
      this.end=xlim[2]
    }
 
    rect(this.start, ypos-height/2, this.end, ypos+height/2, col=this.col, border=this.col)

    if(this.strand=="1"){
      segments(this.start, ypos+height/2, this.start, ypos+height/2+arrowheight)
      arrows(this.start, ypos+height/2+arrowheight, this.start+arrowsize, ypos+height/2+arrowheight, col="gray20", length=0.05, xpd=NA, lwd=1.5)
    } else{
      if(this.strand=="-1"){
        segments(this.end, ypos+height/2, this.end, ypos+height/2+arrowheight)
        arrows(this.end, ypos+height/2+arrowheight, this.end-arrowsize, ypos+height/2+arrowheight, col="gray20", length=0.05, xpd=NA, lwd=1.5)
      } else{
        stop("unknown strand!")
      }
    }

    if(g==focus.gene){
      text(this.name, x=(this.start+this.end)/2, y=ypos+height*0.75, adj=c(0.5, 0), cex=1.1, font=3)
    }
  }

  if(axis){
    axis(side=1, mgp=c(3, 0.5, 0))
  }
}

#####################################################################

plot.annotations.exons<-function(exon.coords, focus.gene, xlim, col.focus="navy", col.other="gray60", axis=T){

  if(!focus.gene%in%exon.coords$id){
    stop(paste("cannot find", focus.gene, "in annotations"))
  }
  ## extract chr
  this.chr=exon.coords$chr[which(exon.coords$geneid==focus.gene)]
  this.exons=exon.coords[which(exon.coords$chr==this.chr & ((exon.coords$start>=xlim[1] & exon.coords$start<=xlim[2]) | (exon.coords$end>=xlim[1] & exon.coords$end<=xlim[2]) | (exon.coords$start<=xlim[1] & exon.coords$end>=xlim[2]))),]

  if(!focus.gene%in%this.genes$id){
   stop(paste(focus.gene, "does not overlap with region")) 
  }

  this.genes$color=col.other
  this.genes$color[which(this.genes$geneid==focus.gene)]=col.focus
  
  print(paste(dim(this.genes)[1], "genes"))
  
  ylim=c(-0.3,0.5)
  ypos=0
  height=0.2

  arrowheight=0.25
  arrowsize=diff(xlim)/30
  
  plot(1, type="n", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, xaxs="i", yaxs="i") ## empty plot
  
  ## abline(h=0, col="gray40")

  for(g in this.genes$id){
    this.strand=this.genes$strand[which(this.genes$geneid==g)]
    this.col=this.genes$color[which(this.genes$id==g)]
    this.start=this.genes$start[which(this.genes$id==g)]
    this.end=this.genes$end[which(this.genes$id==g)]
  
    if(this.start<xlim[1]){
      this.start=xlim[1]
    }

    if(this.end>xlim[2]){
      this.end=xlim[2]
    }

    segments(min(this.start), ypos,  max(this.end), ypos, border=this.col)
    
    rect(this.start, ypos-height/2, this.end, ypos+height/2, col=this.col, border=this.col)
    
    if(this.strand=="+"){
      segments(this.start, ypos+height/2, this.start, ypos+height/2+arrowheight)
      arrows(this.start, ypos+height/2+arrowheight, this.start+arrowsize, ypos+height/2+arrowheight, col="gray20", length=0.075, xpd=NA, lwd=1.5)
    } else{
      segments(this.end, ypos+height/2, this.end, ypos+height/2+arrowheight)
      arrows(this.end, ypos+height/2+arrowheight, this.end-arrowsize, ypos+height/2+arrowheight, col="gray20", length=0.075, xpd=NA, lwd=1.5)
    }
    
  }

  if(axis){
    axis(side=1, mgp=c(3, 0.5, 0))
  }
}

#####################################################################

