###############################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("parameters.R")
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
  
  load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.RData", sep=""))
  load(paste(pathFigures, "RData/data.regland.conservation.RData", sep=""))
  
  regcons=regland.conservation[[sp]][[enh]]
  
  load=FALSE
  
}

#############################################################################################################

if(prepare){
  ## graphical parameters
  distances=c("all", "shortrange", "longrange")
  col.distances=c("black", "steelblue", "indianred")
  names(col.distances)=distances

  smallxdist=c(-0.2, 0, 0.2)
  names(smallxdist)=c("all", "shortrange", "longrange")
    
  prepare=FALSE
}

#############################################################################################################

plot.expdiv.regdiv <- function(regland, feature, expdata, distances, ylab, plot.label, xlab, xax.labels, xax.las){
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
    
    pval = kruskal.test(values, groups)$p.value
    print(paste0("Kruskal-Wallis test : ", ylab, " according to ", feature, " for ", dist, " contacts; p-val:", pval))

    ## compare just first and last class

    firstclass=contact.classes[1]
    lastclass=contact.classes[length(contact.classes)]

    extpval=kruskal.test(values[which(groups%in%c(firstclass, lastclass))], groups[which(groups%in%c(firstclass, lastclass))])$p.value

    print(paste0("Kruskal-Wallis test, first class vs. last class : ", ylab, " according to ", feature, " for ", dist, " contacts; p-val:", extpval))
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
  
  for(dist in distances){
    for(class in contact.classes){
      x=xpos[class]+smallxdist[dist]
      
      med=median.matrix[dist, class]
      ci.low=ci.low.matrix[dist, class]
      ci.high=ci.high.matrix[dist, class]
      
      points(x, med, pch=20, col=col.distances[dist], cex=1.1)
      segments(x, ci.low, x, ci.high,  col=col.distances[dist])
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

selected.distances=c("all", "shortrange", "longrange")

################################################################################################################################

pdf(file=paste(pathFigures, "Figure6.pdf", sep=""), width = 6.85,  height=7)

m=matrix(rep(NA,2*9), nrow=2)
m[1,]=c(rep(c(1,2,3), each=3))
m[2,]=c(rep(4,3), rep(5,2), rep(6,2), rep(7,2))

layout(m)

par(mar = c(6.5, 4.5, 2.5, 0.5))

######################## nb of contacted enhancers and expression characteristics  ########################

## mean expression level as a function of number of contacts
expdata=expdiv[,"MeanRPKM"]
names(expdata)=rownames(expdiv)

plot.expdiv.regdiv(regcons, "class.nb.contacts", expdata, distances=selected.distances, ylab="mean expression level (RPKM)", plot.label="a", xlab="number of contacts", xax.labels=levels(regcons$class.nb.contacts.all), xax.las=2)

###########################################################################
## expression specificity as a function of number of contacts

expdata=expdiv[,paste0("Tau", sp_name)]
names(expdata)=rownames(expdiv)

plot.expdiv.regdiv(regcons, "class.nb.contacts", expdata, distances=selected.distances, ylab="expression specificity", plot.label="b", xlab="number of contacts", xax.labels=levels(regcons$class.nb.contacts.all), xax.las=2)

###########################################################################

## expression conservation as a function of number of contacts

expdata=expdiv[,"CorrectedSpearman"]
names(expdata)=rownames(expdiv)

plot.expdiv.regdiv(regcons, "class.nb.contacts", expdata, distances=selected.distances, ylab="expression conservation", plot.label="c", xlab="number of contacts", xax.labels=levels(regcons$class.nb.contacts.all), xax.las=2)


######################## gene expression profile evolution ################

expdata=expdiv[, "CorrectedSpearman"]
names(expdata)=rownames(expdiv)

## expression conservation as a function of enhancer sequence conservation
plot.expdiv.regdiv(regcons, "class.aln.score", expdata, distances=selected.distances, ylab = "expression conservation", plot.label="d", xlab="enhancer sequence conservation",  xax.labels=levels(regcons$class.aln.score.all), xax.las=2) 


## expression conservation as a function of enhancer synteny conservation
plot.expdiv.regdiv(regcons, "class.synteny.cons", expdata, distances=selected.distances, ylab = "expression conservation", plot.label="e", xlab="synteny conservation", xax.labels=c("<100%", "100%"), xax.las=2)


par(mar = c(6.5, 4.5, 2.5, 0))

## expression conservation as a function of contact conservation
plot.expdiv.regdiv(regcons, "class.contact.cons", expdata, distances=selected.distances, ylab = "expression conservation", plot.label="f", xlab="contact conservation", xax.labels=c("<10%", "10-40%", ">40%"), xax.las=2)

#########################################################################
## empty plot for the legend
par(mar=c(0, 0.1, 0, 0.3)) 
plot.new()

legend("topleft", col=c(col.distances[1], "white", col.distances[2], "white", col.distances[3]), legend = c("all", "", "medium range\n(25 kb - 500 kb)","", "long range\n(> 500 kb)"), box.col="white", bg="white", pch=20, cex=1.1, inset=c(0.1,0.1), xpd=T)

################################################################################################################################

dev.off()

################################################################################################################################
################################################################################################################################
