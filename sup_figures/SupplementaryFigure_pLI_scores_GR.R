###############################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("parameters.R")
}

##############################################################################

if(load){
  
  enh="ENCODE"
  
  load(paste(pathFigures, "RData/data.human.CM2019.AllOrgans.expdiv.RData", sep=""))
  expdiv.human=expdiv
  
  rm("expdiv")
  
  pLI.human = read.table("/home/laverre/Regulatory_landscape/data/pLi_scores/Gene_pLI.txt", h=T)
  pLI.human = pLI.human[!duplicated(pLI.human$GeneID),]
  rownames(pLI.human) = pLI.human$GeneID
  
  common = intersect(rownames(pLI.human), rownames(expdiv.human))

  expdiv.human = expdiv.human[common,]
  expdiv.human$pLI = pLI.human[common,]$pLI
  expdiv.human$class.pLI.all = cut(expdiv.human$pLI, breaks=c(0, 0.1, 0.3, 0.6, 0.9, 1), include.lowest=T, labels=c("<0.1", "0.1-0.3", "0.3-0.6", "0.6-0.9", ">0.9"))
  expdiv.human$gene=expdiv.human$IDHuman
  

  load(paste(pathFigures, "RData/data.regland.conservation.RData", sep=""))
  
  regcons.human=regland.conservation[["human"]][[enh]]

  load=FALSE
  
}

#############################################################################################################

if(prepare){
  ## graphical parameters
  distances=c("all", "shortrange", "longrange")
  col.distances=c("black", "steelblue", "indianred")
  names(col.distances)=distances
  
  smallxdist=c(0, 0.1, 0.2)
  names(smallxdist)=c("all", "shortrange", "longrange")
  
  prepare=FALSE
}

#############################################################################################################

plot.expdiv.regdiv <- function(regland, feature, expdata, distances, ylab, plot.label, smallx.add=0, xlab, xax.labels, xax.las, col=NA, add=FALSE, ylim=NA, mean.plot=F){
  ## go through all classes and distances to compute median and ci
  contact.classes=levels(regland[,paste(feature, distances[1],sep=".")])
  
  median.matrix=matrix(rep(NA, length(distances)*length(contact.classes)), nrow=length(distances))
  rownames(median.matrix)=distances
  colnames(median.matrix)=contact.classes
  
  ci.low.matrix=median.matrix
  ci.high.matrix=median.matrix  
  
  for(dist in distances){
    values <- c()
    groups <- c()
    for(class in contact.classes){
      
      this.genes=regland$gene[which(regland[,paste(feature, dist, sep=".")] == class)]
      b=boxplot(expdata[this.genes], plot=FALSE)
      ci=as.numeric(b$conf)
      
      if (mean.plot == TRUE){
        median.matrix[dist, class]=mean(expdata[this.genes],na.rm=T)
      }else{median.matrix[dist, class]=median(expdata[this.genes],na.rm=T)}

      ci.low.matrix[dist, class]=ci[1]
      ci.high.matrix[dist, class]=ci[2]
      
      values <- c(values, expdata[this.genes])
      groups <- c(groups, rep(class, length(expdata[this.genes])))
    }
  }
  
  ## compute ylim on all values, if not provided
  
  if(any(is.na(ylim))){
    ylim=range(c(as.numeric(ci.low.matrix), as.numeric(ci.high.matrix)))
    addy=diff(ylim)/5
    ylim=ylim+c(-addy, addy)
  }
  
  ## now do the actual plot
  
  xpos=1:length(contact.classes)
  names(xpos)=contact.classes
  xlim=c(0.5, length(contact.classes)+0.5)
  
  cex.mtext = 0.75
  
  if(!add){
    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  } 
  
  for(dist in distances){
    for(class in contact.classes){
      x=xpos[class]+smallxdist[dist]+smallx.add
      
      med=median.matrix[dist, class]
      ci.low=ci.low.matrix[dist, class]
      ci.high=ci.high.matrix[dist, class]
      
      if(!is.na(col)){
        this.col=col
      } else{
        this.col=col.distances[dist]
      }
      points(x, med, pch=20, col=this.col, cex=1.1)
      segments(x, ci.low, x, ci.high,  col=this.col)
    }
  }
  
  if(!add){
    abline(v=xpos[1:length(contact.classes)-1]+diff(xpos)[1]/2, lty=3, col="gray40")
    
    axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1, las=2)
    mtext(ylab, side=2, line=2.95, cex=cex.mtext)
    
    ## plot label
    if(length(xpos)>=5){
      plot.lab.pos=xlim[1]-diff(xlim)/10
    } else{
      plot.lab.pos=xlim[1]-diff(xlim)/10
    }
    
    mtext(plot.label, side=3, line=1.15, at=plot.lab.pos, font=2, cex=1.2)
    
    axis(side=1, cex.axis=1, mgp=c(3, 0.75, 0), at=xpos, labels=rep("",length(xpos)))
    mtext(xax.labels, at=xpos, side=1, line=0.75, cex=0.75, las=xax.las)
    
    mtext(xlab, side=1, line=4.5, cex=cex.mtext)
  }
}

#############################################################################################################
#############################################################################################################

#############################################################################################################

pdf(file=paste(pathFigures, "GenomeResearch_Figures/Figure6_pLI_score.pdf", sep=""), width = 6.85,  height=5)

m=matrix(rep(NA,2*8), nrow=2)
m[1,]=c(1, rep(c(2,3,4), each=2), 5)
m[2,]=c(rep(c(6,7,8,9), each=2))

layout(m)

par(mar = c(6.5, 4.5, 2.5, 0.5))

######################## nb of contacted enhancers and expression characteristics  ########################

#############################################################################

## expression conservation as a function of number of contacts

plot.new()

expdata.human=expdiv.human[,"pLI"]
names(expdata.human)=rownames(expdiv.human)

# expression level
plot.expdiv.regdiv(expdiv.human, "class.MeanRPKM", expdata.human, distances="all", ylab="pLI score",
                   plot.label="A", xlab="mean expression level (RPKM)", xax.labels=levels(expdiv.human$class.MeanRPKM.all),
                   xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(-0.1, 0.8))

# specificity
plot.expdiv.regdiv(expdiv.human, "class.Tau", expdata.human, distances="all", ylab="pLI score",
                   plot.label="B", xlab="expression specifity", xax.labels=levels(expdiv.human$class.Tau.all),
                   xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(-0.1, 0.5))


# nb contact
plot.expdiv.regdiv(regcons.human, "class.nb.contacts", expdata.human, distances="all", ylab="pLI score",
                   plot.label="C", xlab="number of contacts", xax.labels=levels(regcons.human$class.nb.contacts.all),
                   xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0, 0.2))

plot.new()

# sequence
plot.expdiv.regdiv(regcons.human, "class.aln.score", expdata.human, distances="all", ylab="pLI score",
                   plot.label="D", xlab="enhancer sequence conservation", xax.labels=levels(regcons.human$class.aln.score.all),
                   xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0, 0.25))

# phyloP score 
plot.expdiv.regdiv(regcons.human, "class.phyloP.score", expdata.human, distances="all", ylab = "pLI score",
                   plot.label="E", xlab="enhancer phyloP score",  xax.labels=levels(regcons.human$class.phyloP.score.all),
                   xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0, 0.5))

# synteny
plot.expdiv.regdiv(regcons.human, "class.synteny.cons", expdata.human, distances="all", ylab="pLI score",
                   plot.label="F", xlab="synteny conservation rate", xax.labels=levels(regcons.human$class.synteny.cons.all),
                   xax.las=1, smallx.add=-0.15, col="darkred", ylim=c(0, 0.2))

# contact
plot.expdiv.regdiv(regcons.human, "class.contact.cons", expdata.human, distances="all", ylab="pLI score",
                   plot.label="G", xlab="contact conservation rate", xax.labels=levels(regcons.human$class.contact.cons.all),
                   xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0, 0.35))



#############################################################################

dev.off()

#############################################################################
