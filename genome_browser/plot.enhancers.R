#####################################################################

plot.enhancers<-function(enhancer.coords, chr, xlim, col="gray40"){
  nbsamples=length(enhancer.coords)
  ylim=c(0.5, nbsamples+0.5)

  ypos=1:nbsamples
  names(ypos)=names(enhancer.coords)

  height=0.2
  
  plot(1, type="n", xlab="", axes=F, yaxs="i", xaxs="i", ylab="", xlim=xlim, ylim=ylim)
         
  for(sample in names(enhancer.coords)){
    
    this.coords=enhancer.coords[[sample]]
    this.coords=this.coords[which(this.coords$chr==chr & ((this.coords$start>=xlim[1] & this.coords$start<=xlim[2]) | (this.coords$end>=xlim[1] & this.coords$end<=xlim[2]) | (this.coords$start<=xlim[1] & this.coords$end>=xlim[2]))), ]

    nbleft=dim(this.coords)[1]

    print(paste(nbleft, "enhancers in region for ", sample))
    
    this.ypos=ypos[sample]

    rect(this.coords$start, this.ypos-height/2, this.coords$end, this.ypos+height/2, col=col, border=col)
  }
}

#####################################################################
