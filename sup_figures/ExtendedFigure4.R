###############################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

##############################################################################

if(load){
  sp="human"
  
  if(sp == "human"){
    sp_name="Human"
  } else{
    sp_name="Mouse"
  }
  
  enh="ENCODE"
  
  cells <- c("Bcell", "ESC", "adipo")
  col.cells = c("navy", "forestgreen", "darkorange")
  names(col.cells) = cells

  smallxcell=c(-0.2, 0, 0.2)
  names(smallxcell)=cells
  
  load(paste(pathFigures, "RData/data.common.cells.expdiv.RData", sep=""))
  load(paste(pathFigures, "RData/data.common.cells.regland.conservation.RData", sep=""))
  
  regcons=regland.conservation[[sp]][[enh]]
  
  load=FALSE
}

###########################################################################################################

plot.expdiv.regdiv <- function(regland, featurecontact, cells,  expdata, featureexp, ylab, plot.label, xlab, xax.labels, xax.las){
  ## go through all cells to compute median and ci
  firstcell=cells[1]
  contact.classes=levels(regland[[firstcell]][[featurecontact]])
  
  median.matrix=matrix(rep(NA, length(cells)*length(contact.classes)), nrow=length(cells))
  rownames(median.matrix)=cells
  colnames(median.matrix)=contact.classes

  ci.low.matrix=median.matrix
  ci.high.matrix=median.matrix  
  
  for(cell in cells){
    this.regland=regland[[cell]]

    for(class in contact.classes){
      this.genes=this.regland$gene[which(this.regland[, featurecontact] == class)]
      b=boxplot(expdata[this.genes, paste(cell, featureexp, sep="_")], plot=FALSE)
      ci=as.numeric(b$conf)

      median.matrix[cell, class]=median(expdata[this.genes, paste(cell, featureexp, sep="_")],na.rm=T)
      ci.low.matrix[cell, class]=ci[1]
      ci.high.matrix[cell, class]=ci[2]
    }
  }
  
  ## compute ylim on all values
  
  ylim=range(c(as.numeric(ci.low.matrix), as.numeric(ci.high.matrix)))
  addy=diff(ylim)/5
  ylim=ylim+c(-addy, addy)

  ## now do the actual plot
  
  xpos=1:length(contact.classes)
  names(xpos)=contact.classes
  xlim=c(0.5, length(contact.classes)+0.5)
  
  cex.mtext = 0.75
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(cell in cells){
    for(class in contact.classes){
      x=xpos[class]+smallxcell[cell]
      
      med=median.matrix[cell, class]
      ci.low=ci.low.matrix[cell, class]
      ci.high=ci.high.matrix[cell, class]
      
      points(x, med, pch=20, col=col.cells[cell], cex=1.1)
      segments(x, ci.low, x, ci.high,  col=col.cells[cell])
    }
  }
  
  abline(v=xpos[1:length(contact.classes)-1]+diff(xpos)[1]/2, lty=3, col="gray40")
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1, las=2)
  mtext(ylab, side=2, line=2.95, cex=cex.mtext)
  
  ## plot label
  if(length(xpos)==5){
    plot.lab.pos=xlim[1]-diff(xlim)/3.5
  } else{
    plot.lab.pos=xlim[1]-diff(xlim)/2
  }
  
  mtext(plot.label, side=3, line=1.15, at=plot.lab.pos, font=2, cex=1.2)
  
  axis(side=1, cex.axis=1, mgp=c(3, 0.75, 0), at=xpos, labels=rep("",length(xpos)))
  mtext(xax.labels, at=xpos, side=1, line=0.75, cex=0.75, las=xax.las)
  
  mtext(xlab, side=1, line=4.5, cex=cex.mtext)
}

################################################################################################################################
################################################################################################################################

pdf(paste(pathFigures, "/ExtendedFigure4.pdf", sep=""), width=6.85, height=6)

m=matrix(rep(NA,2*9), nrow=2)
m[1,]=c(rep(c(1,2,3), each=3))
m[2,]=c(rep(4,3), rep(5,2), rep(6,2), rep(7,2))

layout(m)

par(mar = c(6.5, 4.5, 2.5, 0.5))

################################################################################################################################

######################## nb of contacted enhancers and expression level  ########################

plot.expdiv.regdiv(regland=regcons, featurecontact="class.nb.contacts", expdata=expdiv_cells, featureexp=paste(sp, "MeanRPKM",sep="_"), cells=cells, ylab="mean expression level (RPKM)", plot.label="a", xlab="number of contacts", xax.labels=levels(regcons[["ESC"]][,"class.nb.contacts"]), xax.las=2)

################################################################################################################################

######################## nb of contacted enhancers and expression conservation, before correction  ########################

plot.expdiv.regdiv(regland=regcons, featurecontact="class.nb.contacts", expdata=expdiv_cells, featureexp="ExpressionConservation", cells=cells, ylab="expression conservation", plot.label="b", xlab="number of contacts", xax.labels=levels(regcons[["ESC"]][,"class.nb.contacts"]), xax.las=2)


######################## nb of contacted enhancers and expression conservation, after correction  ########################

plot.expdiv.regdiv(regland=regcons, featurecontact="class.nb.contacts", expdata=expdiv_cells, featureexp="ResidualExpressionConservation", cells=cells, ylab="exp. cons. (corrected)", plot.label="c", xlab="number of contacts", xax.labels=levels(regcons[["ESC"]][,"class.nb.contacts"]), xax.las=2)

################################################################################################################################

######################## sequence conservation and expression conservation      ########################

plot.expdiv.regdiv(regland=regcons, featurecontact="class.aln.score", expdata=expdiv_cells, featureexp="ResidualExpressionConservation", cells=cells, ylab="exp. cons. (corrected)", plot.label="d", xlab="enhancer sequence conservation", xax.labels=levels(regcons[["ESC"]][,"class.aln.score"]), xax.las=2)

################################################################################################################################

######################## synteny conservation and expression conservation      ########################

plot.expdiv.regdiv(regland=regcons, featurecontact="class.synteny.cons", expdata=expdiv_cells, featureexp="ResidualExpressionConservation", cells=cells, ylab="exp. cons. (corrected)", plot.label="e", xlab="synteny conservation", xax.labels=c("<100%", "100%"), xax.las=2)

################################################################################################################################
#######################  contact conservation and expression conservation        ########################

plot.expdiv.regdiv(regland=regcons, featurecontact="class.contact.cons", expdata=expdiv_cells, featureexp="ResidualExpressionConservation", cells=cells, ylab="exp. cons. (corrected)", plot.label="f", xlab="contact conservation", xax.labels=levels(regcons[["ESC"]][,"class.contact.cons"]), xax.las=2)

################################################################################################################################
# empty plot for the legend
par(mar=c(0, 0.1, 0, 0.3)) 
plot.new()


legend("topleft", bg="white", box.col="white", legend=c("B lymphocytes", "", "embryonic stem cells","", "pre-adipocytes"), col=c(col.cells[1], "white", col.cells[2], "white", col.cells[3]), pch=20, inset=c(0.05, 0.15), xpd=NA)
       
################################################################################################################################

dev.off()

################################################################################################################################

