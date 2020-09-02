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

  for(enh in names(gene.enhancer.contacts[[sp]])){
    real.data=gene.enhancer.contacts[[sp]][[enh]][["real"]]
    nb.contacted.enhancers.real=tapply(real.data$enhancer, as.factor(real.data$gene), function(x) length(unique(x)))
    names(nb.contacted.enhancers.real)=levels(as.factor(real.data$gene))
    
    sim.data=gene.enhancer.contacts[[sp]][[enh]][["simulated"]]
    nb.contacted.enhancers.sim=tapply(sim.data$enhancer, as.factor(sim.data$gene), function(x) length(unique(x)))
    names(nb.contacted.enhancers.sim)=levels(as.factor(sim.data$gene))
    
    nb.contacted.enhancers[[enh]]=list("real"=nb.contacted.enhancers.real, "simulated"=nb.contacted.enhancers.sim)
  }

  ## tissue specificity and other expression statistics

  exp.stats=expstats.cm2019[[sp]]

  ## select genes: protein-coding genes, in the PCHiC data, and in expression stats

  annot=gene.annot[[sp]]
  genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
  genes=intersect(genes, exp.stats$GeneID)

  exp.stats=exp.stats[genes,]

  ## expression specificity
  
  tau=exp.stats$TauLogRPKM
  tau.class=cut(tau, breaks=seq(from=0, to=1, by=0.2), include.lowest=T)
  names(tau.class)=exp.stats$GeneID

  prepare=FALSE
}

##################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##################################################################

pdf(file=paste(pathFigures, "Figure6.pdf", sep=""), width=6.85, height=4)

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

par(mar=c(4.5, 3.5, 2.1, 1.1))
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

    points(x, med, pch=20, col=col.enhancers[enh])
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", 5), cex.axis=0.8)
mtext(c("0-0.2","0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"), at=xpos, side=1, line=0.75, cex=0.8)

mtext("expression specificity index", side=1, line=2.25, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.8)
mtext("number of contacted enhancers", side=2, line=2.5, cex=0.9)

legend("topright", legend=label.enhancers[enhancer.datasets[[sp]]], lty=1, col=col.enhancers[enhancer.datasets[[sp]]], bty="n", cex=0.75, inset=c(0.02, -0.1), seg.len=1,  xpd=NA) 

##################################################################


##################################################################

dev.off()

##################################################################


