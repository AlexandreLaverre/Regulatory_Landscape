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

    obs$nb_celltypes <- apply(obs[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
    obs$celltype_class<- cut(obs$nb_celltypes, breaks=c(0:7, max(obs$nb_celltypes)), include.lowest=T)
    levels(obs$celltype_class)=c(as.character(1:7), ">7")

    sim$nb_celltypes <- apply(sim[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
    sim$celltype_class<- cut(sim$nb_celltypes, breaks=c(0:7, max(sim$nb_celltypes)), include.lowest=T)
    levels(sim$celltype_class)=c(as.character(1:7), ">7")

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

pdf(paste(pathFigures, "SupplementaryFigure10.pdf", sep=""), width=4.49, height=3)

m=matrix(c(1,2), nrow=1)

layout(m)

##########################################################################

for(sp in c("human", "mouse")){
  obs=all.data[[sp]][["obs"]]
  sim=all.data[[sp]][["sim"]]

  classes=levels(obs$celltype_class)
  nbclass=length(classes)
  
  xpos=1:nbclass
  smallx=c(-0.25, 0.25)

  ylim=c(0, max(c(obs$NbSamplesExpressed, sim$NbSamplesExpressed)))

  plot(1, type="n", xlab="", ylab="", axes=F)
  for(i in 1:nbclass){
    wobs=which(obs$celltype_class==classes[i])
    boxplot(obs$NbSamplesExpressed[wobs], outline=F, add=T, axes=F, at=xpos[i]+smallx[1], border=dataset.colors["Original"])
    
    wsim=which(sim$celltype_class==classes[i])
    boxplot(sim$NbSamplesExpressed[wsim], outline=F, add=T, axes=F, at=xpos[i]+smallx[1], border=dataset.colors["Original"])
    
  }
  
}

##########################################################################

dev.off()

##########################################################################
