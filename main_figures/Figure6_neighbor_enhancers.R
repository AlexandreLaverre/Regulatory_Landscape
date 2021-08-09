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

  load(paste(pathFigures, "RData/data.mouse.CM2019.AllOrgans.expdiv.RData", sep=""))
  expdiv.mouse=expdiv

  rm("expdiv")
  
  load(paste(pathFigures, "RData/data.regland.conservation.neighbor.enhancers.RData", sep=""))
  
  regcons.human=regland.conservation[["human"]][[enh]]
  regcons.mouse=regland.conservation[["mouse"]][[enh]]
  
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

plot.expdiv.regdiv <- function(regland, feature, expdata, distances, ylab, plot.label, smallx.add=0, xlab, xax.labels, xax.las, col=NA, add=FALSE, ylim=NA){
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

      median.matrix[dist, class]=median(expdata[this.genes],na.rm=T)
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
      plot.lab.pos=xlim[1]-diff(xlim)/3.5
    } else{
      plot.lab.pos=xlim[1]-diff(xlim)/2
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

pdf(file=paste(pathFigures, "GenomeResearch_Figures/Figure6_neighbor_enhancers.pdf", sep=""), width = 6.85,  height=7)

m=matrix(rep(NA,2*9), nrow=2)
m[1,]=c(rep(c(1,2,3), each=3))
m[2,]=c(rep(4,3), rep(5,2), rep(6,4))

layout(m)

par(mar = c(6.5, 4.5, 2.5, 0.5))

######################## nb of contacted enhancers and expression characteristics  ########################

## mean expression level as a function of number of contacts
expdata.mouse=expdiv.mouse[,"MeanRPKM"]
names(expdata.mouse)=rownames(expdiv.mouse)

expdata.human=expdiv.human[,"MeanRPKM"]
names(expdata.human)=rownames(expdiv.human)

plot.expdiv.regdiv(regcons.human, "class.nb.contacts", expdata.human, distances="all",  ylab="mean expression level (RPKM)", plot.label="A", xlab="number of contacts", xax.labels=levels(regcons.human$class.nb.contacts.all), xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(7, 9.5))

plot.expdiv.regdiv(regcons.mouse, "class.nb.contacts", expdata.mouse, distances="all", add=T, smallx.add=0.15, col="navy")

## tests 

common.human=intersect(rownames(regcons.human), names(expdata.human))
kw.human=kruskal.test(expdata.human[common.human]~regcons.human[common.human, "class.nb.contacts.all"])
print(paste("human, class nb contacts vs. mean RPKM, pval", kw.human$p.value))


common.mouse=intersect(rownames(regcons.mouse), names(expdata.mouse))
kw.mouse=kruskal.test(expdata.mouse[common.mouse]~regcons.mouse[common.mouse, "class.nb.contacts.all"])
print(paste("mouse, class nb contacts vs. mean RPKM, pval", kw.mouse$p.value))

###########################################################################
## expression specificity as a function of number of contacts

expdata.mouse=expdiv.mouse[,"TauMouse"]
names(expdata.mouse)=rownames(expdiv.mouse)

expdata.human=expdiv.human[,"TauHuman"]
names(expdata.human)=rownames(expdiv.human)

plot.expdiv.regdiv(regcons.human, "class.nb.contacts", expdata.human, distances="all", ylab="expression specificity", plot.label="B", xlab="number of contacts", xax.labels=levels(regcons.human$class.nb.contacts.all), xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0.65, 0.78))
plot.expdiv.regdiv(regcons.mouse, "class.nb.contacts", expdata.mouse, distances="all",  add=T, smallx.add=0.15, col="navy")


## tests 

common.human=intersect(rownames(regcons.human), names(expdata.human))
kw.human=kruskal.test(expdata.human[common.human]~regcons.human[common.human, "class.nb.contacts.all"])
print(paste("human, class nb contacts vs. Tau, pval", kw.human$p.value))


common.mouse=intersect(rownames(regcons.mouse), names(expdata.mouse))
kw.mouse=kruskal.test(expdata.mouse[common.mouse]~regcons.mouse[common.mouse, "class.nb.contacts.all"])
print(paste("mouse, class nb contacts vs. Tau, pval", kw.mouse$p.value))

#############################################################################

## expression conservation as a function of number of contacts

expdata.mouse=expdiv.mouse[,"CorrectedSpearman"]
names(expdata.mouse)=rownames(expdiv.mouse)

expdata.human=expdiv.human[,"CorrectedSpearman"]
names(expdata.human)=rownames(expdiv.human)

plot.expdiv.regdiv(regcons.human, "class.nb.contacts", expdata.human, distances="all", ylab="expression conservation", plot.label="C", xlab="number of contacts", xax.labels=levels(regcons.human$class.nb.contacts.all), xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(-0.05, 0.1))
plot.expdiv.regdiv(regcons.mouse, "class.nb.contacts", expdata.mouse, distances="all", add=T, smallx.add=0.15, col="navy")


## tests 

common.human=intersect(rownames(regcons.human), names(expdata.human))
kw.human=kruskal.test(expdata.human[common.human]~regcons.human[common.human, "class.nb.contacts.all"])
print(paste("human, class nb contacts vs. expression conservation, pval", kw.human$p.value))


common.mouse=intersect(rownames(regcons.mouse), names(expdata.mouse))
kw.mouse=kruskal.test(expdata.mouse[common.mouse]~regcons.mouse[common.mouse, "class.nb.contacts.all"])
print(paste("mouse, class nb contacts vs. expression conservation, pval", kw.mouse$p.value))


########################## gene expression profile evolution ################

## expression conservation as a function of enhancer sequence conservation
plot.expdiv.regdiv(regcons.human, "class.aln.score", expdata.human, distances="all", ylab = "expression conservation", plot.label="D", xlab="enhancer sequence conservation",  xax.labels=levels(regcons.human$class.aln.score.all), xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(-0.06, 0.12))
plot.expdiv.regdiv(regcons.mouse, "class.aln.score", expdata.mouse, distances="all",add=T,smallx.add=0.15, col="navy")

## tests 

common.human=intersect(rownames(regcons.human), names(expdata.human))
kw.human=kruskal.test(expdata.human[common.human]~regcons.human[common.human, "class.aln.score.all"])
print(paste("human, class aln score vs. expression conservation, pval", kw.human$p.value))


common.mouse=intersect(rownames(regcons.mouse), names(expdata.mouse))
kw.mouse=kruskal.test(expdata.mouse[common.mouse]~regcons.mouse[common.mouse, "class.aln.score.all"])
print(paste("mouse, class aln score vs. expression conservation, pval", kw.mouse$p.value))

#########################################################################

## ##  expression conservation as a function of enhancer synteny conservation
##  plot.expdiv.regdiv(regcons.human, "class.synteny.cons", expdata.human, distances="all", ylab = "expression conservation", plot.label="E", xlab="synteny conservation", xax.labels=c("<100%", "100%"), xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(0.00, 0.06))
## plot.expdiv.regdiv(regcons.mouse, "class.synteny.cons", expdata.mouse, distances="all", add=T, smallx.add=0.15, col="navy")


## ## tests 

## common.human=intersect(rownames(regcons.human), names(expdata.human))
## kw.human=kruskal.test(expdata.human[common.human]~regcons.human[common.human, "class.synteny.cons.all"])
## print(paste("human, class synteny vs. expression conservation, pval", kw.human$p.value))


## common.mouse=intersect(rownames(regcons.mouse), names(expdata.mouse))
## kw.mouse=kruskal.test(expdata.mouse[common.mouse]~regcons.mouse[common.mouse, "class.synteny.cons.all"])
## print(paste("mouse, class synteny vs. expression conservation, pval", kw.mouse$p.value))

#########################################################################
## expression conservation as a function of contact conservation


#############################################################################

## empty plot for the legend

##legend("bottomright", col=c("darkred", "white", "navy"), legend = c("human", "", "mouse"), box.col="white", bg="white", pch=20, cex=1.1, inset=c(-0.4,0.1), xpd=T)

#############################################################################

dev.off()

#############################################################################
#############################################################################
