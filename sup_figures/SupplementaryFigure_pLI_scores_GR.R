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
  
  load(paste(pathFigures, "RData/data.regland.conservation.phyloP.RData", sep=""))
  
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

pdf(file=paste(pathFigures, "GenomeResearch_Figures/Figure6_pLI_vs_phyloP.pdf", sep=""), width = 5.85,  height=6)

m=matrix(rep(NA,2*6), nrow=2)
m[1,]=c(rep(c(1,2), each=3))
m[2,]=c(rep(c(3,4), each=3))

layout(m)

par(mar = c(6.5, 4.5, 2.5, 0.5))

######################## nb of contacted enhancers and expression characteristics  ########################

#############################################################################

## expression conservation as a function of number of contacts

# nb contact
# expdata.human=regcons.human[,"nb.contacts.all"]
# names(expdata.human)=rownames(regcons.human)
# 
# plot.expdiv.regdiv(expdiv.human, "class.pLI", expdata.human, distances="all", ylab="number of contacts",
#                    plot.label="A", xlab="pLI scores", xax.labels=levels(expdiv.human$class.pLI.all),
#                    xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(25, 40))
# 
# # sequence
# expdata.human=regcons.human[,"mean.aln.score.all"]
# names(expdata.human)=rownames(regcons.human)
# 
# plot.expdiv.regdiv(expdiv.human, "class.pLI", expdata.human, distances="all", ylab="enhancer sequence conservation",
#                    plot.label="B", xlab="pLI scores", xax.labels=levels(expdiv.human$class.pLI.all),
#                    xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0.44, 0.5))
# 
# # synteny
# expdata.human=regcons.human[,"fr.synteny.cons.all"]
# names(expdata.human)=rownames(regcons.human)
# 
# plot.expdiv.regdiv(expdiv.human, "class.pLI", expdata.human, distances="all", ylab="synteny conservation rate",
#                    plot.label="C", xlab="pLI scores", xax.labels=levels(expdiv.human$class.pLI.all), mean.plot=TRUE,
#                    xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0.9, 1))
# 
# # contact
# expdata.human=regcons.human[,"fr.contact.cons.all"]
# names(expdata.human)=rownames(regcons.human)
# 
# plot.expdiv.regdiv(expdiv.human, "class.pLI", expdata.human, distances="all", ylab="contact conservation rate",
#                    plot.label="D", xlab="pLI scores", xax.labels=levels(expdiv.human$class.pLI.all),
#                    xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0.14,0.3))
# 
# 
# # sequence
# expdata.human=regcons.human[,"mean.phyloP.score.all"]
# names(expdata.human)=rownames(regcons.human)
# 
# plot.expdiv.regdiv(expdiv.human, "class.pLI", expdata.human, distances="all", ylab="enhancer phyloP score",
#                    plot.label="B", xlab="pLI scores", xax.labels=levels(expdiv.human$class.pLI.all),
#                    xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0.44, 0.5))

# class pLI vs phyloP score
expdata.human=regcons.human[,"mean.aln.score.all"]
names(expdata.human)=rownames(regcons.human)

plot.expdiv.regdiv(expdiv.human, "class.pLI", expdata.human, distances="all", ylab="enhancer phyloP score",
                   plot.label="B", xlab="pLI scores", xax.labels=levels(expdiv.human$class.pLI.all),
                   xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0.44, 0.5))


# class phyloP score vs pLI vs 
expdata.human=expdiv.human[,"pLI"]
names(expdata.human)=rownames(expdiv.human)

plot.expdiv.regdiv(regcons.human, "class.phyloP.score", expdata.human, distances="all", ylab = "pLI score",
                   plot.label="B", xlab="enhancer phyloP score",  xax.labels=levels(regcons.human$class.phyloP.score.all),
                   xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0, 0.5))


#############################################################################

dev.off()

#############################################################################
