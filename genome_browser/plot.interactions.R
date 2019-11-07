######################################################################

plot.interactions<-function(interactions, xlim, focus.bait, col.contact="indianred", col.otherbait="red"){

  this.int=interactions[which(interactions$bait_ID%in%focus.bait),]
  
  ylim=c(-1, 1)

  plot(1,type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

  height=0.5
  ypos=0
  
  if(dim(this.int)[1]>0){
    col=rep(col.contact, dim(this.int)[1])
    col[which(this.int$baited_frag=="baited")]=col.otherbait
    rect(this.int$start, ypos-height/2, this.int$end, ypos+height/2, col=col, border=NA)
  }

}

######################################################################
