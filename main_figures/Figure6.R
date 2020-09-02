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
  tau.class=cut(tau, breaks=seq(from=0, to=1, by=0.25), include.lowest=T)
  names(tau.class)=exp.stats$GeneID

  prepare=FALSE
}

##################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##################################################################

pdf(file=paste(pathFigures, "Figure6.pdf", sep=""), width=4.49, height=4)

##################################################################

## classes of expression specificity

xpos=c(1,2,3,4)
names(xpos)=levels(tau.class)

smallx=c(-0.3, -0.1, 0.1, 0.3)
names(smallx)=enh.datasets[[sp]]

ylim=c(0, 50)
xlim=c(0, 5)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

for(enh in enh.datasets[[sp]]){
  nb.contacts=nb.contacted.enhancers[[sp]][["real"]]
    
  for(class in levels(tau.class)){
    this.genes=names(tau.class)[which(tau.class==class)]
    this.genes=intersect(this.genes, names(nb.contacts))
    
    x=xpos[class]+smallx[enh]
        
    boxplot(nb.contacts[this.genes], at=x, notch=T, outline=F, add=T, axes=F)
  }
}

##################################################################

dev.off()

##################################################################


