#######################################################################################

options(stringsAsFactors = FALSE)

source("parameters.R")

#######################################################################################

print("loading data")

load(paste(pathFigures, "RData/data.neighbor.enhancers.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
load(paste(pathFigures, "RData/data.bait.annotation.RData", sep=""))
load(paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))

print("done")

## we load sequence and synteny conservation later

#######################################################################################

## minDistance and maxDistance defined in parameters.R

## min and max number of enhancers for contacts and synteny conservation

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

  for(enh in enhancer.datasets[[ref]]){
  
    print(paste(ref, enh))
    
    ## we start from the gene-enhancer contact data
    
    contacts=gene.enhancer.contacts[[ref]][[enh]][["real"]]
      
    ## we select only previously filtered orthologous genes
    
    contacts=contacts[which(contacts$gene%in%ortho[,ref]),]

    ## we add sequence conservation column

    load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref,"2", tg,".RData", sep=""))
    contacts$align_score=pcungapped[contacts$enhancer]
   
    
    ## prepare result list

    results=list("gene"=all.genes)
    
    ## nb enhancers by gene
    
    nb.contacts=as.numeric(table(factor(contacts$gene, levels=all.genes)))
    
    ## divide nb of contacts into 5 classes
    
    results[["nb.contacts.all"]]=nb.contacts
    results[["class.nb.contacts.all"]]=cut(nb.contacts, breaks=c(1, 10, 20, 30, 40, max(nb.contacts)), include.lowest=T, labels=c("1-10", "11-20", "21-30", "30-40", ">40"))
    
    ## mean conservation score by gene
    
    mean.aln.score=tapply(contacts$align_score, factor(contacts$gene, levels=all.genes), mean, na.rm=T)
    
    results[["mean.aln.score.all"]]=mean.aln.score
    results[["class.aln.score.all"]]=cut(mean.aln.score, breaks=seq(from=0, to=1, length=6), include.lowest=T)
    
    
    results=as.data.frame(results)
    
    regland.conservation[[ref]][[enh]]=results
  }
}

#######################################################################################

save(regland.conservation, file=paste(pathFigures, "RData/data.regland.conservation.neighbor.enhancers.RData",sep=""))

#######################################################################################
#######################################################################################
