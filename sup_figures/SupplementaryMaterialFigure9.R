#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

##########################################################################

if(load){
 
  load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.bait.annotation.RData", sep=""))
  load(paste(pathFigures, "RData/data.gene.expression.RData", sep=""))

  minRPKM=1

  maxval=c(7, 5)
  names(maxval)=c("human", "mouse")

  labels=c("a", "b")
  names(labels)=c("human", "mouse")
  
  load=FALSE
}

##########################################################################

if(prepare){
  ## determine number of cell types in which contacts were observed

  all.data=list()

  for(sp in c("human","mouse")){
    obs=observed.contacts[[sp]]
    sim=simulated.contacts[[sp]]
    
    info=sampleinfo[[sp]]
    rownames(info)=info$Sample.ID
    
    samples=info$Sample.ID 
    celltypes=info$Broad.cell.type.or.tissue
    names(celltypes)=samples

    M=maxval[sp]
    
    obs$nb_celltypes <- apply(obs[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
    sim$nb_celltypes <- apply(sim[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
    
    ## bait annotation

    bait.annot=bait.info[[sp]]
    obs$geneID=bait.annot[obs$id_bait, "gene_ID"]
    sim$geneID=bait.annot[sim$id_bait, "gene_ID"]
    
    ## expression

    exp=avgexp.cm2019[[sp]]
   
    nbsamples.exp=apply(exp, 1, function(x) length(which(x>=minRPKM)))
    names(nbsamples.exp)=rownames(exp)

    ## select baits associated with a single gene, in expression data

    obs=obs[which(obs$geneID%in%rownames(exp)),]
    sim=sim[which(sim$geneID%in%rownames(exp)),]

    ## we compute the maximum number of samples in which a gene has chromatin contacts

    nb.celltypes.bygene.obs=tapply(obs$nb_celltypes, as.factor(obs$geneID), max)
    nb.celltypes.bygene.sim=tapply(sim$nb_celltypes, as.factor(sim$geneID), max)

    fac.celltypes.bygene.obs=cut(nb.celltypes.bygene.obs, breaks=c(0:M, max(obs$nb_celltypes)), include.lowest=T)
    levels(fac.celltypes.bygene.obs)=c(as.character(1:M), paste0(">",M))

    fac.celltypes.bygene.sim=cut(nb.celltypes.bygene.sim, breaks=c(0:M, max(sim$nb_celltypes)), include.lowest=T)
    levels(fac.celltypes.bygene.sim)=c(as.character(1:M), paste0(">",M))

    res.obs=data.frame("gene"=names(nb.celltypes.bygene.obs), "nbcontact"= nb.celltypes.bygene.obs, "classcontact"=fac.celltypes.bygene.obs, "nbexp"=nbsamples.exp[names(nb.celltypes.bygene.obs)])
    res.sim=data.frame("gene"=names(nb.celltypes.bygene.sim), "nbcontact"= nb.celltypes.bygene.sim, "classcontact"=fac.celltypes.bygene.sim, "nbexp"=nbsamples.exp[names(nb.celltypes.bygene.sim)])
   
    all.data[[sp]]=list("obs"=res.obs, "sim"=res.sim)
  }
  
  prepare=FALSE
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##########################################################################

pdf(paste(pathFigures, "SupplementaryMaterialFigure9.pdf", sep=""), width=6.85, height=3.5)

m=matrix(c(rep(1, 7), rep(2, 7)), nrow=1)

layout(m)

##########################################################################

for(sp in c("human", "mouse")){
  exp=avgexp.cm2019[[sp]]
  nbsamples.tot=dim(exp)[2]
  
  obs=all.data[[sp]][["obs"]]
  sim=all.data[[sp]][["sim"]]

  classes=levels(obs$classcontact)
  nbclass=length(classes)
  
  xpos=1:nbclass
  smallx=c(-0.25, 0.25)

  xlim=c(0.5, 9.5)
  ylim=c(0, 100)

  par(mar=c(3.5, 3.5, 1.5, 0.5))
      
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)
  
  for(i in 1:nbclass){
    wobs=which(obs$classcontact==classes[i])
    boxplot(100*obs$nbexp[wobs]/nbsamples.tot, add=T, axes=F, at=xpos[i]+smallx[1], col=dataset.colors["Original"], notch=T, pch=20)
    
    wsim=which(sim$classcontact==classes[i])
    boxplot(100*sim$nbexp[wobs]/nbsamples.tot, add=T, axes=F, at=xpos[i]+smallx[2], col=dataset.colors["Simulated"], notch=T, pch=20)
    
  }

  axis(side=1, at=1:nbclass, labels=classes, cex.axis=0.9, mgp=c(3, 0.5, 0))
  mtext("number of cell types w. chromatin contacts", side=1, line=2, cex=0.7, at=(nbclass+1)/2)
  
  axis(side=2, cex.axis=0.9, mgp=c(3, 0.75, 0))
  mtext("percentage of samples w. detectable expression", side=2, line=2, cex=0.7)

  mtext(labels[sp], side=3, at=xlim[1]-diff(xlim)/6.5, line=0, font=2)
  
  if (sp == "mouse"){
    mtext(sp, side=3, cex=0.75, at=3.5)
    legend("topright", legend=c("PCHi-C data", "simulated data"), fill=dataset.colors[c("Original", "Simulated")], bty='n', inset=c(-0.02, -0.01), xpd=NA)
  } else{
    mtext(sp, side=3, cex=0.75, at=4.5)
  }
}

##########################################################################

dev.off()

##########################################################################
