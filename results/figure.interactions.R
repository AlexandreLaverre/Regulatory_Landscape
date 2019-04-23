###################################################################

library(DescTools)

###################################################################

plot.arc <- function(x1, x2, y, col, lwd, plot=TRUE){
  x0=(x1+x2)/2
  R=abs(x2-x1)/sqrt(2)
  y0=y-abs(x2-x1)/2
  
  ymax=y0+R

  if(plot==TRUE){
    DrawArc(x=x0, y=y0, rx=R, theta.1=pi/4, theta.2=3*pi/4, col=col, lwd=lwd)
  }
  else{
    return(list("x0"=x0, "y0"=y0, "ymax"=ymax, "R"=R))
  }
}

###################################################################

plot.interactions<-function(gene, interactions, coords, x.scale="kb"){
  this.inter=interactions[which(interactions$Ensembl.Gene.ID==gene),]
  this.inter=this.inter[which(this.inter$chr==this.inter$chr_bait), ]## only intra-chromosome
  
  this.coords=coords[which(coords$id==gene),]
  this.strand=this.coords$strand
  this.chr=this.coords$chr

  print(paste(dim(this.inter)[1], "interactions"))
  
  xlim=range(c(this.coords$start, this.coords$end, this.inter$start_bait, this.inter$end_bait, this.inter$start, this.inter$end))
  print(xlim)

  all.coords=coords[which(coords$chr==this.coords$chr),]
  intersect.regions=unlist(lapply(1:dim(all.coords)[1], function(i) {M=max(xlim[1], all.coords$start[i]); m=min(xlim[2], all.coords$end[i]); return(M<=m)}))
  all.coords=all.coords[which(intersect.regions),]

  print(paste(dim(all.coords)[1], "genes"))
  
  all.genes=all.coords$id
  col.genes=rep("gray60", length(all.genes))
  names(col.genes)=all.genes
  col.genes[gene]="red"
  
  layout(matrix(c(1, 1, 2), nrow=3))


  ## first plot: interactions

  this.inter$x1=(this.inter$start_bait+this.inter$end_bait)/2
  this.inter$x2=(this.inter$start+this.inter$end)/2

  ## estimate y values to see the arcs

  ypos=0.5

  ymax=max(unlist(lapply(1:dim(this.inter)[1], function(i) plot.arc(x1=this.inter$x1[i], x2=this.inter$x2[i], y=ypos, plot=FALSE)[["ymax"]])))

  ylim=c(0, ymax)

  par(mar=c(0.5,2,0.5,1))
  plot(1, type="n", xlab="", ylab="",axes=F, xlim=xlim, ylim=ylim)

  for(i in 1:dim(this.inter)[1]){
    plot.arc(x1=this.inter$x1[i], x2=this.inter$x2[i], y=ypos, col="red", lwd=1)
  }


  abline(h=ypos, lty=1, col="gray40", lwd=0.5)
  
  
  ## second plot: annotations
  
  par(mar=c(1.5,2,0.15,1))
  plot(1, type="n", xlab="", ylab="",axes=F, xlim=xlim, ylim=c(0,1))
  ypos=c(0.35,0.75)
  
  if(this.strand=="-1"){
    names(ypos)=c("1", "-1")
  } else{
    names(ypos)=c("-1", "1")
  }
  
  height=0.075
  
  abline(h=mean(c(ypos[2], ypos[1]+height)), lty=2, col="gray80")
  
  arrow.width=diff(xlim)/100
  arrow.height=0.075
  
  for(g in all.genes){
    this.strand=as.character(all.coords$strand[which(all.coords$id==g)])
    this.ypos=ypos[this.strand]
    this.start=all.coords$start[which(all.coords$id==g)]
    this.end=all.coords$end[which(all.coords$id==g)]
    
    rect(this.start, this.ypos-height, this.end, this.ypos+height, col=col.genes[g], border=col.genes[g]) 
    
    if(this.strand=="1"){
      segments(this.start, this.ypos+height, this.start, this.ypos+height+arrow.height, col=col.genes[g])
      arrows(this.start, this.ypos+height+arrow.height, this.start+arrow.width, this.ypos+height+arrow.height, length=0.05, col=col.genes[g])
    } else{
      segments(this.end, this.ypos+height, this.end, this.ypos+height+arrow.height, col=col.genes[g])
      arrows(this.end, this.ypos+height+arrow.height, this.end-arrow.width, this.ypos+height+arrow.height, length=0.05, col=col.genes[g])
    }

    if(g==gene){
      text(this.coords$name, x=(this.start+this.end)/2, y=this.ypos++height+3.5*arrow.height, xpd=NA, cex=0.8)
    }
  }

  if(x.scale=="kb"){
    xax=pretty(xlim/1000)
    axis(side=1, at=xax*1000, labels=paste(xax, "kb"), mgp=c(3,0.5,0))
    mtext(paste("chr", this.chr, sep=""), at=xlim[1], line=0.5, cex=0.7, side=1)
  } else{
    if(x.scale=="Mb"){
      xax=pretty(xlim/1000)
      axis(side=1, at=xax*1000, labels=paste(xax, "Mb"), mgp=c(3,0.5,0))
      mtext(paste("chr", this.chr, sep=""), at=xlim[1], line=0.5, cex=0.7, side=1)
    }  else{
      axis(side=1, mgp=c(3,0.5,0))
      mtext(paste("chr", this.chr, sep=""), at=xlim[1], line=0.5, cex=0.7, side=1)
    }
  }
}

###################################################################

## load=TRUE

if(load==TRUE){
  
  pathInteractions="/pandata/alaverre/data/human_promoter-other_hg38_annot5kb_filtred.txt"
  inter=read.table(pathInteractions, sep="\t", h=T, stringsAsFactors=F)
  
  pathAnnot="/pandata/necsulea/RegulatoryLandscapes/data/gene_annotations/human/"
  
  gene.coords=read.table(paste(pathAnnot,"GeneInfo_Ensembl91.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  
  colnames(gene.coords)=c("id", "biotype", "description", "chr", "start", "end", "strand")
  rownames(gene.coords)=gene.coords$id
  
  gene.names=read.table(paste(pathAnnot,"GeneNames_Ensembl91.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(gene.names)=gene.names$Gene.stable.ID

  gene.coords$name=gene.names[gene.coords$id, "Gene.name"]

  load=FALSE
}

###################################################################

## plot.interactions("ENSG00000164690", inter, gene.coords)


plot.interactions(inter$Ensembl.Gene.ID[1], inter, gene.coords)

###################################################################

