##################################################################

source("parameters.R")

##################################################################

if(load){

  sp="human"

  load(paste(pathFigures, "RData/data.gene.expression.RData", sep=""))
  load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
  load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
  
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

  prepare=FALSE
}

##################################################################



##################################################################




##################################################################


