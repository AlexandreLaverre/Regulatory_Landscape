######################################################################

plot.interactions<-function(interactions, xlim, focus.bait, col.contact="indianred", col.otherbait="red"){

  ylim=c(-1, 1)

  plot(1,type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

  height=0.5
  ypos=0
  
  if(dim(interactions)[1]>0){
    col=rep(col.contact, dim(interactions)[1])
    col[which(interactions$baited_frag=="baited")]=col.otherbait
    rect(interactions$start, ypos-height/2, interactions$end, ypos+height/2, col=col, border=NA)
  }

}

######################################################################
