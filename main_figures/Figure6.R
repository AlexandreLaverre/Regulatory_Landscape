##################################################################

source("parameters.R")

##################################################################

if(load){

  sp="human"

  load(paste(pathFigures, "RData/data.gene.expression.RData", sep=""))
  load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
  load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
  
  load=FALSE
}

##################################################################

if(prepare){

  ## compute number of contacted enhancers per gene

  nb.contacted.enhancers=list()
  median.distance=list()

  for(enh in names(gene.enhancer.contacts[[sp]])){
    real.data=gene.enhancer.contacts[[sp]][[enh]][["real"]]
    nb.contacted.enhancers.real=tapply(real.data$enhancer, as.factor(real.data$gene), function(x) length(unique(x)))
    names(nb.contacted.enhancers.real)=levels(as.factor(real.data$gene))

    real.mediandist=tapply(real.data$dist, as.factor(real.data$gene), median)
    names(real.mediandist)=levels(as.factor(real.data$gene))

    real.data.multi=real.data[which(real.data$nb_sample>=2),]
    nb.contacted.enhancers.real.multi=tapply(real.data.multi$enhancer, as.factor(real.data.multi$gene), function(x) length(unique(x)))
    names(nb.contacted.enhancers.real.multi)=levels(as.factor(real.data.multi$gene))
    
    ## simulated data
    
    sim.data=gene.enhancer.contacts[[sp]][[enh]][["simulated"]]
    nb.contacted.enhancers.sim=tapply(sim.data$enhancer, as.factor(sim.data$gene), function(x) length(unique(x)))
    names(nb.contacted.enhancers.sim)=levels(as.factor(sim.data$gene))

    sim.mediandist=tapply(sim.data$dist, as.factor(sim.data$gene), median)
    names(sim.mediandist)=levels(as.factor(sim.data$gene))
    
     sim.data.multi=sim.data[which(sim.data$nb_sample>=2),]
    nb.contacted.enhancers.sim.multi=tapply(sim.data.multi$enhancer, as.factor(sim.data.multi$gene), function(x) length(unique(x)))
    names(nb.contacted.enhancers.sim.multi)=levels(as.factor(sim.data.multi$gene))

    ## results
    
    nb.contacted.enhancers[[enh]]=list("real"=nb.contacted.enhancers.real, "simulated"=nb.contacted.enhancers.sim, "real.multi"=nb.contacted.enhancers.real.multi, "sim.multi"=nb.contacted.enhancers.sim.multi)
    median.distance[[enh]]=list("real"=real.mediandist, "simulated"=sim.mediandist)
  }

  ## tissue specificity and other expression statistics

  exp.stats=expstats.cm2019[[sp]]
  avgexp.tissage=avgexp.cm2019[[sp]]

  ## select genes expressed in at least one sample
  exp.stats=exp.stats[which(exp.stats$MaxRPKM>0),]

  ## select genes: protein-coding genes, in the PCHiC data, and in expression stats

  annot=gene.annot[[sp]]
  genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
  genes=intersect(genes, exp.stats$GeneID)

  exp.stats=exp.stats[genes,]
  avgexp.tissage=avgexp.tissage[genes,]

  ## expression specificity classes
  
  tau=exp.stats$TauLogRPKM
  tau.class=cut(tau, breaks=seq(from=0, to=1, by=0.2), include.lowest=T)
  names(tau.class)=exp.stats$GeneID

  ## average expression class

  avgexp=log2(exp.stats$AverageRPKM+1)
  avgexp.class=cut(avgexp, breaks=quantile(avgexp, p=seq(from=0, to=1, by=0.2)), include.lowest=T)
  names(avgexp.class)=exp.stats$GeneID

  ## age order

  if(sp=="human"){
    exp.stats$MaxAge[which(exp.stats$MaxAge%in%c("oldTeenager", "youngTeenager"))]="teenager"
    exp.stats$MaxAge[which(exp.stats$MaxAge%in%c("Senior"))]="senior"
    exp.stats$MaxAge[which(exp.stats$MaxAge%in%c("youngAdult"))]="adult"
    exp.stats$MaxAge[which(exp.stats$MaxAge%in%c("youngMidage", "oldMidage"))]="midage"

    age.avgexp=unlist(lapply(colnames(avgexp.tissage), function(x) unlist(strsplit(x, split="\\."))[2]))
    age.avgexp[which(age.avgexp%in%c("oldTeenager", "youngTeenager"))]="teenager"
    age.avgexp[which(age.avgexp%in%c("Senior"))]="senior"
    age.avgexp[which(age.avgexp%in%c("youngAdult"))]="adult"
    age.avgexp[which(age.avgexp%in%c("youngMidage", "oldMidage"))]="midage"
    
    age.order=c(paste(4:20, "wpc", sep=""), "newborn", "infant", "toddler", "school", "teenager", "adult", "midage", "senior")
  }

  if(sp=="mouse"){
    age.order=c(paste("e", 10:18, sep=""), paste("P", 0:63, sep=""))
  }

  age.order=intersect(age.order, exp.stats$MaxAge)
    
  prepare=FALSE
}

##################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##################################################################

pdf(file=paste(pathFigures, "Figure6.pdf", sep=""), width=6.85, height=5.5)

m=matrix(rep(NA, 2*10), nrow=2)
m[1,]=c(rep(1, 5), rep(2,5))
m[2,]=c(rep(3, 10))

layout(m)

##################################################################

## classes of expression specificity

xpos=c(1,2,3,4,5)
names(xpos)=levels(tau.class)

smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]

ylim=c(0, 50)
xlim=c(0.5, 5.5)

par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  nb.contacts=nb.contacted.enhancers[[enh]][["real"]]
  
  for(class in levels(tau.class)){
    this.genes=names(tau.class)[which(tau.class==class)]
    this.genes=intersect(this.genes, names(nb.contacts))
    
    x=xpos[class]+smallx[enh]
    
    b=boxplot(nb.contacts[this.genes], plot=FALSE)
    med=median(nb.contacts[this.genes])
    ci=as.numeric(b$conf)

    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:4]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", 5), cex.axis=0.8)
mtext(c("0-0.2","0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"), at=xpos, side=1, line=1, cex=0.8)

mtext("expression specificity index", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("nb. contacted enhancers", side=2, line=2.5, cex=0.9)

legend("topright", legend=label.enhancers[enhancer.datasets[[sp]]], lty=1, col=col.enhancers[enhancer.datasets[[sp]]], bty="o", cex=1.1, inset=c(-0.01, -0.25), seg.len=1,  xpd=NA, box.col="white", bg="white")

mtext("A", side=3, at=-0.2, font=2, cex=1.1, line=1.5)

##################################################################

## average expression level

xpos=c(1,2,3,4,5)
names(xpos)=levels(avgexp.class)

smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]

ylim=c(0, 50)
xlim=c(0.5, 5.5)

par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  nb.contacts=nb.contacted.enhancers[[enh]][["real"]]
  
  for(class in levels(avgexp.class)){
    this.genes=names(avgexp.class)[which(avgexp.class==class)]
    this.genes=intersect(this.genes, names(nb.contacts))
    
    x=xpos[class]+smallx[enh]
    
    b=boxplot(nb.contacts[this.genes], plot=FALSE)
    med=median(nb.contacts[this.genes])
    ci=as.numeric(b$conf)

    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:4]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", 5), cex.axis=0.8)
mtext(c("low", "medium", "high"), at=xpos[c(1,3,5)], side=1, line=1, cex=0.8)
mtext("expression level class", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("nb. contacted enhancers", side=2, line=2.5, cex=0.9)

mtext("B", side=3, at=-0.2, font=2, cex=1.1, line=1.5)

##################################################################

## age of maximum expression vs. number of enhancers

xpos=1:length(age.order)
names(xpos)=age.order

smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]

ylim=c(0, 50)
xlim=c(0.5, length(age.order)+0.5)

par(mar=c(6.1, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

minRPKM=5

for(enh in enhancer.datasets[[sp]]){
  nb.contacts=nb.contacted.enhancers[[enh]][["real"]]
  
  for(class in age.order){
    this.genes=exp.stats$GeneID[which(exp.stats$MaxAge==class & exp.stats$MaxRPKM>1)]
    
    x=xpos[class]+smallx[enh]
    
    b=boxplot(nb.contacts[this.genes], plot=FALSE)
    med=median(nb.contacts[this.genes], na.rm=T)
    ci=as.numeric(b$conf)

    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:(length(age.order)-1)]+0.5, lty=3, col="gray40")

abline(v=xpos[14]+0.5, lty=2)

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(age.order)), cex.axis=0.8)
mtext(age.order, at=xpos, side=1, cex=0.75, las=2, line=1)
mtext("developmental stage with maximum expression", side=1, line=4.9, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("nb. contacted enhancers", side=2, line=2.5, cex=0.9)

mtext("C", side=3, at=-0.85, font=2, cex=1.1, line=1.5)

##################################################################

dev.off()

##################################################################


