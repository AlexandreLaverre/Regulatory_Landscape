#####################################################################

plot.annotations.genes<-function(gene.coords, focus.gene, gene.biotypes="all", xlim, col.focus="navy", col.other="gray60", axis=T, axisunit=NA, cex.name=1){

  if(!focus.gene%in%gene.coords$id){
    stop(paste("cannot find", focus.gene, "in annotations"))
  }
  ## extract chr
  this.chr=gene.coords$chr[which(gene.coords$id==focus.gene)]
  this.genes=gene.coords[which(gene.coords$chr==this.chr & ((gene.coords$start>=xlim[1] & gene.coords$start<=xlim[2]) | (gene.coords$end>=xlim[1] & gene.coords$end<=xlim[2]) | (gene.coords$start<=xlim[1] & gene.coords$end>=xlim[2]))),]

  if(!all(gene.biotypes=="all")){
    this.genes=this.genes[which(this.genes$biotype%in%gene.biotypes),]
  }
  
  if(!focus.gene%in%this.genes$id){
   stop(paste(focus.gene, "does not overlap with region")) 
  }

  this.genes$color=col.other
  this.genes$color[which(this.genes$id==focus.gene)]=col.focus
  
  print(paste(dim(this.genes)[1], "genes"))
  
  ylim=c(-3,3)
  ypos.bystrand=c(1,-1)
  names(ypos.bystrand)=c("+", "-")
  height=0.65

  arrowheight=0.65
  arrowsize=diff(xlim)/75
  arrowlength=0.035
  
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
 
    rect(this.start, ypos-height/2, this.end, ypos+height/2, col=this.col, border=this.col, xpd=NA)

    if(this.strand=="+"){
      segments(this.start, ypos+height/2, this.start, ypos+height/2+arrowheight, col=this.col)
      arrows(this.start, ypos+height/2+arrowheight, this.start+arrowsize, ypos+height/2+arrowheight, length=arrowlength, xpd=NA, lwd=1.5, col=this.col)
      if(g==focus.gene){
        text(this.name, x=(this.start+this.end)/2, y=ypos+height*0.75, adj=c(0.5, 0), cex=cex.name, font=3)
      }
    } else{
      if(this.strand=="-"){
        segments(this.end, ypos+height/2, this.end, ypos+height/2+arrowheight, col=this.col)
        arrows(this.end, ypos+height/2+arrowheight, this.end-arrowsize, ypos+height/2+arrowheight, length=arrowlength, xpd=NA, lwd=1.5, col=this.col)

         if(g==focus.gene){
           text(this.name, x=(this.start+this.end)/2, y=ypos-height*0.75, adj=c(0.5, 1), cex=cex.name, font=3)
         }
        
      } else{
        print(this.strand)
        stop("unknown strand!")
      }
    }

   
  }

  if(axis){
    if(is.na(axisunit)){
       axis(side=1, mgp=c(3, 0.5, 0))
     } else{
       if(axisunit%in%c("kb", "Kb", "k", "K")){
         xax=pretty(xlim/1e3)
         xval=paste(xax, "kb", sep="")
         axis(side=1, at=xax*1e3, labels=xval, mgp=c(3, 0.5, 0))
       } else{
         if(axisunit%in%c("Mb", "mb", "m", "M")){
           xax=pretty(xlim/1e6)
           xval=paste(xax, "Mb", sep="")
           axis(side=1, at=xax*1e6, labels=xval, mgp=c(3, 0.5, 0))
         }
       }
     }
  }
}

#####################################################################
