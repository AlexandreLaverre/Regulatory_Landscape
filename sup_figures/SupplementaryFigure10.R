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
    obs$celltype_class<- cut(obs$nb_celltypes, breaks=c(0:M, max(obs$nb_celltypes)), include.lowest=T)
    levels(obs$celltype_class)=c(as.character(1:M), paste0(">",M))

    sim$nb_celltypes <- apply(sim[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
    sim$celltype_class<- cut(sim$nb_celltypes, breaks=c(0:M, max(sim$nb_celltypes)), include.lowest=T)
    levels(sim$celltype_class)=c(as.character(1:M), paste0(">",M))

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

    obs$NbSamplesExpressed=nbsamples.exp[obs$geneID]
    sim$NbSamplesExpressed=nbsamples.exp[sim$geneID]

    all.data[[sp]]=list("obs"=obs, "sim"=sim)
  }
  
  prepare=FALSE
}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##########################################################################

pdf(paste(pathFigures, "SupplementaryFigure10.pdf", sep=""), width=6.85, height=3.5)

m=matrix(c(rep(1, 7), rep(2, 6)), nrow=1)

layout(m)

##########################################################################

for(sp in c("human", "mouse")){
  obs=all.data[[sp]][["obs"]]
  sim=all.data[[sp]][["sim"]]

  classes=levels(obs$celltype_class)
  nbclass=length(classes)
  
  xpos=1:nbclass
  smallx=c(-0.25, 0.25)

  xlim=c(0.5, nbclass+0.5)
  
  if (sp == "human"){ylim=c(70,120)}else{ylim=c(50,100)}

  par(mar=c(3.5, 3.5, 1.5, 1.1))
      
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)
  
  for(i in 1:nbclass){
    wobs=which(obs$celltype_class==classes[i])
    obs_stats = boxplot(obs$NbSamplesExpressed[wobs], plot=F)
    
    wsim=which(sim$celltype_class==classes[i])
    sim_stats= boxplot(sim$NbSamplesExpressed[wobs], plot=F)
    
    
    points(i, obs_stats$stats[3,], pch=20, col=dataset.colors["Original"])
    segments(i, obs_stats$conf[1,], i, obs_stats$conf[2,], col=dataset.colors["Original"])
    
    points(i, sim_stats$stats[3,], pch=20, col=dataset.colors["Simulated"])
    segments(i, sim_stats$conf[1,], i, sim_stats$conf[2,], col=dataset.colors["Simulated"])
  }

  mtext(sp, side=3, cex=0.75)

  axis(side=1, at=1:nbclass, labels=classes, cex.axis=0.9, mgp=c(3, 0.5, 0))
  mtext("number of cell types in which contacts are observed", side=1, line=2, cex=0.7)
  
  axis(side=2, cex.axis=0.9, mgp=c(3, 0.75, 0))
  mtext("number of samples with detectable expression", side=2, line=2, cex=0.7)

  mtext(labels[sp], side=3, at=xlim[1]-diff(xlim)/(nbclass-0.8), line=0, font=2)
  
  if (sp == "human"){legend("topleft", legend=c("PCHi-C data", "simulated data"), col=dataset.colors[c("Original", "Simulated")], pch=20,
                            bty='n', inset=c(0.01, 0.05), xpd=NA)}
}

##########################################################################

dev.off()

##########################################################################
