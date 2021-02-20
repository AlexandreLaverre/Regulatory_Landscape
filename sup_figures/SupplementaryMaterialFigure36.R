######################################################################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  library(Hmisc)
  source("../main_figures/parameters.R")
}

##############################################################################

if(load){
  sp="human"
  
  if (sp == "human"){
    sp_name="Human"
  } else{
    sp_name="Mouse"
  }
  
  load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.RData", sep=""))
  expdiv$EuclideanSimilarity = 1-expdiv$EuclideanDistance

  load(paste(pathFigures, "RData/data.regland.conservation.RData", sep=""))
  
  regcons=regland.conservation[[sp]]

  load=FALSE
}


#############################################################################################################

if(prepare){
  ## graphical parameters
  enhancers=enhancer.datasets[[sp]]
  
  smallxenh=c(-0.15, -0.075, 0.075, 0.15)
  names(smallxenh)=enhancers
    
  prepare=FALSE
}

################################################################################################################################
############################## Cardoso-Moreira  - Euclidean Similarity ##########################################################

plot.expdiv.regdiv <- function(regland, feature, expdata, distance, enhancer.list, ylab, plot.label, xlab, xax.labels, xax.las){
  ## go through all classes and enhancer datasets to compute median and ci

  enh1=enhancer.list[1]
  contact.classes=levels(regland[[enh1]][,paste(feature, distance, sep=".")])
  
  median.matrix=matrix(rep(NA, length(enhancer.list)*length(contact.classes)), nrow=length(enhancer.list))
  rownames(median.matrix)=enhancer.list
  colnames(median.matrix)=contact.classes

  ci.low.matrix=median.matrix
  ci.high.matrix=median.matrix  
  
  for(enh in enhancer.list){
    for(class in contact.classes){
      this.regland=regland[[enh]]
      this.genes=this.regland$gene[which(this.regland[,paste(feature, distance, sep=".")] == class)]
      b=boxplot(expdata[this.genes], plot=FALSE)
      ci=as.numeric(b$conf)

      median.matrix[enh, class]=median(expdata[this.genes],na.rm=T)
      ci.low.matrix[enh, class]=ci[1]
      ci.high.matrix[enh, class]=ci[2]
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
  
  for(enh in enhancers){
    for(class in contact.classes){
      x=xpos[class]+smallxenh[enh]
      
      med=median.matrix[enh, class]
      ci.low=ci.low.matrix[enh, class]
      ci.high=ci.high.matrix[enh, class]
      
      points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, ci.low, x, ci.high,  col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:length(contact.classes)-1]+diff(xpos)[1]/2, lty=3, col="gray40")
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1, las=2)
  mtext(ylab, side=2, line=3.2, cex=cex.mtext)
  
  ## plot label
  if(length(xpos)==5){
    plot.lab.pos=xlim[1]-diff(xlim)/5.5
  }

  if(length(xpos)==2){
    plot.lab.pos=xlim[1]-diff(xlim)/3.8
  }

  if(length(xpos)==3){
    plot.lab.pos=xlim[1]-diff(xlim)/5.1
  }

  mtext(plot.label, side=3, line=0.5, at=plot.lab.pos, font=2, cex=1.2)
  
  axis(side=1, cex.axis=1, mgp=c(3, 0.75, 0), at=xpos, labels=rep("",length(xpos)))
  mtext(xax.labels, at=xpos, side=1, line=0.75, cex=0.75, las=xax.las)

  if(xax.las==2){
    mtext(xlab, side=1, line=4.5, cex=cex.mtext)
  }

  if(xax.las==1){
    mtext(xlab, side=1, line=2.5, cex=cex.mtext)
  }
}

################################################################################################################################

pdf(file=paste(pathFigures, "SupplementaryMaterialFigure36.pdf", sep=""), width=6.85, height=11)
m=matrix(rep(NA, 4*8), nrow=4)

m[1,]=c(rep(1, 4), rep(5, 4))
m[2,]=c(rep(2, 4), rep(6, 4))
m[3,]=c(rep(3, 3), rep(9, 1), rep(7, 3), rep(10, 1))
m[4,]=c(rep(4, 4), rep(8, 4))


layout(m)


################################################################################################################################

par(mar=c(4.75, 5.25, 2.5, 1.5)) # bottom, left, top, right

## Euclidean similarity, uncorrected

expdata=expdiv[,"EuclideanSimilarity"]
names(expdata)=rownames(expdiv)

## number of contacts
plot.expdiv.regdiv(regcons, "class.nb.contacts", expdata, "all", enhancers, ylab="1-Euclidean distance", plot.label="a", xlab="number of contacts", xax.labels=levels(regcons[[enhancers[1]]]$class.nb.contacts.all), xax.las=1)

legend("bottomright", legend=enhancers, pch=20, col=col.enhancers, cex=1, bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))

## enhancer sequence conservation

plot.expdiv.regdiv(regcons, "class.aln.score", expdata, "all", enhancers, ylab="1-Euclidean distance", plot.label="b", xlab="enhancer sequence conservation", xax.labels=levels(regcons[[enhancers[1]]]$class.aln.score.all), xax.las=1)

## synteny conservation

plot.expdiv.regdiv(regcons, "class.synteny.cons", expdata, "all", enhancers, ylab="1-Euclidean distance", plot.label="c", xlab="synteny conservation", xax.labels=c("<100%", "100%"), xax.las=1)

## contact conservation
par(mar=c(4.75, 5.25, 2.5, 2.5)) # bottom, left, top, right

plot.expdiv.regdiv(regcons, "class.contact.cons", expdata, "all", enhancers, ylab="1-Euclidean distance", plot.label="d", xlab="contact conservation", xax.labels=c("<10%", "10-40%", ">40%"), xax.las=1)

################################################################################################################################

par(mar=c(4.75, 5.25, 2.5, 1.5)) # bottom, left, top, right

## Euclidean similarity, corrected

expdata=expdiv[,"CorrectedEuclideanSimilarity"]
names(expdata)=rownames(expdiv)

## number of contacts
plot.expdiv.regdiv(regcons, "class.nb.contacts", expdata, "all", enhancers, ylab="1-Euclidean distance\n(corrected)", plot.label="e", xlab="number of contacts", xax.labels=levels(regcons[[enhancers[1]]]$class.nb.contacts.all), xax.las=1)

## enhancer sequence conservation

plot.expdiv.regdiv(regcons, "class.aln.score", expdata, "all", enhancers, ylab="1-Euclidean distance\n(corrected)", plot.label="f", xlab="enhancer sequence conservation", xax.labels=levels(regcons[[enhancers[1]]]$class.aln.score.all), xax.las=1)

## synteny conservation

plot.expdiv.regdiv(regcons, "class.synteny.cons", expdata, "all", enhancers, ylab="1-Euclidean distance\n(corrected)", plot.label="g", xlab="synteny conservation", xax.labels=c("<100%", "100%"), xax.las=1)


## contact conservation

par(mar=c(4.75, 5.25, 2.5, 2.5)) # bottom, left, top, right

plot.expdiv.regdiv(regcons, "class.contact.cons", expdata, "all", enhancers, ylab="1-Euclidean distance\n(corrected)", plot.label="h", xlab="contact conservation", xax.labels=c("<10%", "10-40%", ">40%"), xax.las=1)

################################################################################################################################

dev.off()

################################################################################################################################
