##########################################################################

options(stringsAsFactors = FALSE)

#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
  
  set.seed(19)
}

##################################################################

if(load){
  load(paste(pathFigures, "RData/data.bootstrap.chicago.scores.RData",sep=""))
  
  load=FALSE
}

#####################################################################################
## make the actual figure

for(sp in c("human", "mouse")){
   
  if (sp == "human"){
    pdf.name = "GenomeResearch_Figures/Supplemental_Fig_S9.pdf"
  }else{
    pdf.name = "GenomeResearch_Figures/SupplementaryMaterialFigure21.pdf"
  }
  
  tg=setdiff(c("human", "mouse"), sp)
  sampleinfo.tg=sampleinfo[[tg]]
  
  ###################################################################################
  
  pdf(paste(pathFigures, pdf.name, sep=""), width=6.85, height=5.5)
  
  par(mai = c(0.5, 0.5, 0.5, 0.2)) # bottom, left, top, right
  mtext.CEX = 0.75
  
  m=matrix(rep(NA, 2*10), nrow=2)
  m[1,]=c(rep(1,5), rep(2,5))
  m[2,]=c(rep(3,5), rep(4,5))
  layout(m)
  
############################## CHICAGO score according to dist ###############################
  
  if(sp == "human"){
    ylim=c(6, 10)
  }else{
    ylim=c(6, 11)
  }
  
  plot(1, type="n", xlim=c(0.5, 40.5), ylim=ylim, xlab="", ylab="", axes=F)
  xpos=1:39
  
  for(enh in enhancer.datasets[[sp]]){
    lines(xpos, chicago.dist[[sp]][[enh]][xpos], col=col.enhancers[enh])
    segments(xpos, chicago.dist.conf.low[[sp]][[enh]][xpos], xpos, chicago.dist.conf.high[[sp]][[enh]][xpos], col=col.enhancers[enh])
  }
  
  xax=c(1, 10,  20, 30, 40)
  axlab=as.character(c(0.05,  0.5,  1, 1.5, 2))
  axis(side=1, mgp=c(3, 0.65, 0), at=xax, labels=axlab)
  mtext("distance to promoters (Mb)", side=1, line=2, cex=mtext.CEX)
  
  axis(side=2, mgp=c(3, 0.75, 0))
  mtext("CHiCAGO score", side=2, line=2.5, cex=mtext.CEX)
  
  legend("topleft", legend=label.enhancers[enhancer.datasets[[sp]]], lty=1, 
         col=col.enhancers[enhancer.datasets[[sp]]], bty="n", cex=1, seg.len=1)

  ## plot label
  mtext("A", side=3, font=2, at=-7.2, line=2.5)
  
  
  ############################## Sequence conservation ###############################
  xpos=seq(1, nb_chicago_class, 1)
  names(xpos) = 1:nb_chicago_class
  smallx=c(-0.15, -0.075, 0.075, 0.15)
  names(smallx)=enhancer.datasets[[sp]]
  
  if (sp == "human"){
    ylim=c(25, 55)
  }else{
    ylim=c(35, 55)
  }
  
  xlim=c(0.5, nb_chicago_class+0.5)
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(row in 1:length(cons.seq[[sp]][[enh]])){
      
      x=xpos[row]+smallx[enh]
      points(x, cons.seq[[sp]][[enh]][row], pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, cons.seq.conf.low[[sp]][[enh]][row], x, cons.seq.conf.high[[sp]][[enh]][row], col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:nb_chicago_class-1]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", nb_chicago_class))
  mtext(1:nb_chicago_class, at=xpos, side=1, line=1, cex=mtext.CEX)
  mtext("CHiCAGO score decile", side=1, line=2.5, cex=mtext.CEX)
  
  axis(side=2, mgp=c(3, 0.75, 0))
  mtext(paste0("% aligned sequence in ", tg), side=2, line=2.5, cex=mtext.CEX)

  ## plot label
  mtext("B", side=3, font=2, at=-0.9, line=2.5)
  
  ############################## Synteny conservation ###############################
  if (sp == "human"){
    ylim=c(93, 100)
  }else{
    ylim=c(92, 100)
  }
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(row in 1:length(cons.synteny[[sp]][[enh]])){
      
      x=xpos[row]+smallx[enh]
      points(x, cons.synteny[[sp]][[enh]][row], pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, cons.synteny.conf.low[[sp]][[enh]][row], x, cons.synteny.conf.high[[sp]][[enh]][row], col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:nb_chicago_class-1]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", nb_chicago_class))
  mtext(1:nb_chicago_class, at=xpos, side=1, line=1, cex=mtext.CEX)
  mtext("CHiCAGO score decile", side=1, line=2.5, cex=mtext.CEX)
  
  axis(side=2, mgp=c(3, 0.75, 0))
  mtext(paste0("% pairs in conserved synteny in ", tg), side=2, line=2.5, cex=mtext.CEX)
  
  ## plot label
  mtext("C", side=3, font=2, at=-0.9, line=2.5)
  
############################## Contact conservation ###############################
  
  if (sp == "human"){
    ylim=c(10, 60)
  }else{
    ylim=c(20, 55)
  }
    
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(row in 1:length(cons.contact[[sp]][[enh]])){
      
      x=xpos[row]+smallx[enh]
      points(x, cons.contact[[sp]][[enh]][row], pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, cons.contact.conf.low[[sp]][[enh]][row], x, cons.contact.conf.high[[sp]][[enh]][row], col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:nb_chicago_class-1]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", nb_chicago_class))
  mtext(1:nb_chicago_class, at=xpos, side=1, line=1, cex=mtext.CEX)
  mtext("CHiCAGO score decile", side=1, line=2.5, cex=mtext.CEX)
  
  axis(side=2, mgp=c(3, 0.75, 0))
  mtext(paste0("% pairs in conserved contact in ", tg), side=2, line=2.5, cex=mtext.CEX)

  
  ## plot label
  mtext("D", side=3, font=2, at=-0.9, line=2.5)
  
  #####################################################################################
  
  dev.off()
  
  #####################################################################################
  
}

