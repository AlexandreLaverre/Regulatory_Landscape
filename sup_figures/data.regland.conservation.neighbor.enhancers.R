#######################################################################################

options(stringsAsFactors = FALSE)

source("../main_figures/parameters.R")

#######################################################################################

print("loading data")

load(paste(pathFigures, "RData/data.neighbor.enhancers.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
load(paste(pathFigures, "RData/data.bait.annotation.RData", sep=""))
load(paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))

print("done")

## we load sequence and synteny conservation later

#######################################################################################

dist.classes=c("all",  "neighbors")
min.distances=c(minDistance, minDistance)
max.distances=c(maxDistance, 100000)

## min and max number of enhancers for synteny conservation

min.nb=5
max.nb=100

#######################################################################################

regland.conservation=list()

for(ref in c("human", "mouse")){
  
  tg=setdiff(c("human", "mouse"), ref)

  regland.conservation[[ref]]=list()

  ## full list of genes: filtered ortho genes, baited, baits have filtered contacts
  
  all.genes=ortho[,ref]

  frag.contacts=observed.contacts[[ref]]
  filtered.baits=unique(frag.contacts$id_bait)
  this.bait.info=bait.info[[ref]]
  this.bait.info=this.bait.info[which(this.bait.info$ID%in%filtered.baits),]
  baited.genes=unique(unlist(lapply(this.bait.info$gene_ID, function(x) unlist(strsplit(x, split=",")))))
  
  all.genes=intersect(all.genes, baited.genes)

  enh="ENCODE"
  
  print(paste(ref, enh))
  
  ## we start from the gene-enhancer data
  
  contacts=gene.enhancer.contacts[[ref]][[enh]][["real"]]
  contacts$id=paste(contacts$GeneID, contacts$EnhancerID, sep="-")
    
  ## we select only previously filtered orthologous genes
  
  contacts=contacts[which(contacts$GeneID%in%ortho[,ref]),]
  
  ## same minimum distance as for the PCHi-C contacts
  
  contacts=contacts[which(contacts$Distance>=minDistance),]
  
  ## sequence conservation column
  
  load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref,"2", tg,".RData", sep=""))
  contacts$align_score=pcungapped[contacts$EnhancerID]

  ## synteny conservation
  load(paste(pathFigures, "RData/data.synteny.conservation.neighbor.enhancers.RData",sep=""))
  synteny=synteny.cons[[paste(ref, "2", tg, sep="")]]
  synteny$id=paste(synteny$GeneID, synteny$EnhancerID, sep="-")
  
  ## prepare result list
  
  results=list("gene"=all.genes)
  
  for(i in 1:length(dist.classes)){
    
    ## select interactions in a given distance range
    
    dist.class=dist.classes[i]
    min.dist=min.distances[i]
    max.dist=max.distances[i]
    
    print(dist.class)
    
    filtered.contacts=contacts[which(contacts$Distance>=min.dist & contacts$Distance<=max.dist),]
    
    ## mean conservation score by gene
    
    mean.aln.score=tapply(filtered.contacts$align_score, factor(filtered.contacts$GeneID, levels=all.genes), mean, na.rm=T)
    
    results[[paste("mean.aln.score", dist.class, sep=".")]]=mean.aln.score
    results[[paste("class.aln.score", dist.class, sep=".")]]=cut(mean.aln.score, breaks=seq(from=0, to=1, length=6), include.lowest=T)
    
    ## synteny conservation
    
    filtered.synteny=synteny[which(synteny$id%in%filtered.contacts$id),]
    
    fr.cons.synt=tapply(filtered.synteny$ConservedSynteny,  factor(filtered.synteny$GeneID, levels=all.genes), function(x) length(which(x=="yes"))/length(x))
    
    nb.enh.synt=as.numeric(table(factor(filtered.synteny$GeneID, levels=all.genes)))
    
    ## if fewer than min.nb enhancers, we assign NA values to synteny conservation
    
    fr.cons.synt[which(nb.enh.synt < min.nb | nb.enh.synt > max.nb)]=NA
    
    results[[paste("fr.synteny.cons",dist.class,sep=".")]]=fr.cons.synt
    
    results[[paste("class.synteny.cons", dist.class, sep=".")]]=factor(fr.cons.synt<1,levels=c(TRUE, FALSE), labels=c("broken", "conserved"))
    
  }
  
  results=as.data.frame(results)
  
  regland.conservation[[ref]]=results
  
}

#######################################################################################

save(regland.conservation, file=paste(pathFigures, "RData/data.regland.conservation.neighbor.enhancers.RData",sep=""))

#######################################################################################
#######################################################################################
