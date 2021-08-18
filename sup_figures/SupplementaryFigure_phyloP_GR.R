#########################################################################################################################

objects=ls()

if(!"pathFigures"%in%objects){
  
  source("parameters.R")
  
  library(ape)
  library(vioplot)
  
  set.seed(19)
  
  load=T
}

#########################################################################################################################

if(load){
  ref_sp = "human"
  target_sp = "mouse"
  
  enhancers = enhancer.datasets[[ref_sp]]
  
  selenh="ENCODE"
  load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", ref_sp, ".RData", sep=""))
  
  
  ## expression divergence
  load(paste(pathFigures, "RData/data.human.CM2019.AllOrgans.expdiv.RData", sep=""))
  expdiv.human=expdiv
  
  load(paste(pathFigures, "RData/data.mouse.CM2019.AllOrgans.expdiv.RData", sep=""))
  expdiv.mouse=expdiv
  
  rm("expdiv")
  
  load(paste(pathFigures, "RData/data.regland.conservation.RData", sep=""))
  
  regcons.human=regland.conservation[["human"]][[selenh]]
  regcons.mouse=regland.conservation[["mouse"]][[selenh]]
  
  distances=c("all")
  col.distances=c("black", "steelblue", "indianred")
  names(col.distances)=distances
  
  smallxdist=c(0, 0.1, 0.2)
  names(smallxdist)=c("all")
  
  load=F
}


#########################################################################################################################

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
    mtext(xax.labels, at=xpos, side=1, line=0.75, cex=0.7, las=xax.las)
    
    mtext(xlab, side=1, line=4.5, cex=cex.mtext)
  }
}

#########################################################################################################################
## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#########################################################################################################################

pdf(paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_phastCons.pdf", sep=""), width=6, height=6)

m=matrix(rep(NA,6*4), nrow=4)
m[1,]=c(rep(c(1,2), each=2), 5, 5)
m[2,]=c(rep(c(1,2), each=2), 6, 6)
m[3,]=c(rep(c(3,4), each=2), 6, 6)
m[4,]=c(rep(c(3,4), each=2), 7, 7)

layout(m)

######################## phyloP scores vs distance to promoters  ########################

par(mai = c(0.6, 0.2, 0.1, 0.2)) #bottom, left, top and right
par(mar=c(4.1, 4, 2.5, 1.5))

nbclasses=length(levels(frag_align_obs$dist_class))
xpos=1:nbclasses

xlim=c(-0.5, max(xpos)+1)

## axis position
class_leg <- c("0", "0.5", "1", "1.5", "2")
xax=seq(from=0, to=max(xpos)+1, by=10)

n=1

for (ref_sp in c("human", "mouse")){
  for(type in c("restriction fragments", "enhancers")){
    load(paste(pathFigures, "RData/data.bootstrap.conservation.distance.",type,".",ref_sp,".phastCons_score.RData", sep=""))
    
    ylim=range(c(ci.low.obs, ci.high.obs, ci.low.sim, ci.high.sim))
    
    dy=diff(ylim)/20
    ylim=ylim+c(-dy, dy)
    
    plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")
    
    points(xpos, mean.val.obs, col=dataset.colors["Original"], pch=20)
    segments(xpos, ci.low.obs, xpos, ci.high.obs, col=dataset.colors["Original"])
    
    points(xpos, mean.val.sim, col=dataset.colors["Simulated"], pch=20)
    segments(xpos, ci.low.sim, xpos, ci.high.sim, col=dataset.colors["Simulated"])
    
    axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)
    
    if(type=="enhancers"){
      mtext("distance to promoters (Mb)", side=1, line=2.2, cex=0.8)
    }
    
    if(type=="restriction fragments"){
      mtext("distance to baits (Mb)", side=1, line=2.2, cex=0.8)
    }
    
    axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
    mtext("phastCons score", side=2, line=3, cex=0.8)
    
    mtext(paste(ref_sp,", ", type, sep=""), side=3, cex=0.7, line=1)
    
    mtext(LETTERS[n], side=3, line=1, at=-7.75, font=2, cex=1.1)
    n = n+1
  }

}


#######################################################################################################
par(mai = c(0.5, 0.2, 0.2, 0.5)) #bottom, left, top and right

plot(1, type="n", xlab="", ylab="", axes=F, main="", xaxs="i")
legend("left", fill=dataset.colors, border=dataset.colors, legend = c("PCHi-C data", "simulated data"), bty='n', cex=1.3, xpd=T, inset=c(-0.01, -0.1), horiz=FALSE)

par(mai = c(0.5, 0.7, 0.2, 0.1)) #bottom, left, top and right

########################## gene expression profile evolution ################
expdata.mouse=expdiv.mouse[,"CorrectedSpearman"]
names(expdata.mouse)=rownames(expdiv.mouse)

expdata.human=expdiv.human[,"CorrectedSpearman"]
names(expdata.human)=rownames(expdiv.human)

## expression conservation as a function of enhancer sequence conservation

plot.expdiv.regdiv(regcons.human, "class.phastCons.score", expdata.human, distances="all", ylab = "expression conservation",
                   plot.label="E", xlab="phastCons scores",  xax.labels=levels(regcons.human$class.phastCons.score.all),
                   xax.las=2, smallx.add=-0.15, col="darkred", ylim=c(-0.03, 0.2))

plot.expdiv.regdiv(regcons.mouse, "class.phastCons.score", expdata.mouse, distances="all",add=T,smallx.add=0.15, col="navy")

legend("topleft", col=c("darkred", "white", "navy"), legend = c("human", "", "mouse"), box.col="white", bg="white", pch=20, cex=1,
       inset=c(0.05,0))

#######################################################################################################
dev.off()

#######################################################################################################
