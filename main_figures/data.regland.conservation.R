#######################################################################################

options(stringsAsFactors = FALSE)

source("parameters.R")

#######################################################################################

print("loading data")

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData", sep=""))
load(paste(pathFigures, "RData/data.bait.annotation.RData", sep=""))
load(paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))

print("done")

## we load sequence and synteny conservation later

#######################################################################################

## minDistance and maxDistance defined in parameters.R

dist.classes=c("all", "shortrange", "longrange")
min.distances=c(minDistance, minDistance, 500001)
max.distances=c(maxDistance, 500000, maxDistance)

## min and max number of enhancers for contacts and synteny conservation

min.nb=5
max.nb=100

#######################################################################################

regland.conservation=list()

for(ref in c("human", "mouse")){
  
  tg=setdiff(c("human", "mouse"), ref)

  sampleinfo.ref=sampleinfo[[ref]]
  sampleinfo.tg=sampleinfo[[tg]]

  samples.ref=sampleinfo.ref[,"Sample.ID"]
  samples.tg=sampleinfo.tg[,"Sample.ID"]

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

    ## we load synteny conservation - already filtered, min sequence conservation >= 10% threshold
    
    load(paste(pathFigures, "/RData/data.synteny.conservation.", ref, ".RData", sep=""))
    
    synteny=conserv_synteny[[enh]][[tg]][["synt_obs"]]
    
    ## we load contact conservation - already filtered, min sequence conservation 0.4

    contact.cons=contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]][["obs"]]

    contact.cons$is_conserved=apply(contact.cons[,samples.tg], 1, function(x) any(x>0))
    
    ## prepare result list

    results=list("gene"=all.genes)
    
    for(i in 1:length(dist.classes)){
      
      ## select interactions in a given distance range
      
      dist.class=dist.classes[i]
      min.dist=min.distances[i]
      max.dist=max.distances[i]
      
      print(dist.class)

      filtered.contacts=contacts[which(contacts$dist>=min.dist & contacts$dist<=max.dist),]

      ## nb enhancers by gene

      nb.contacts=as.numeric(table(factor(filtered.contacts$gene, levels=all.genes)))

      ## divide nb of contacts into 5 classes
     
      results[[paste("nb.contacts",dist.class,sep=".")]]=nb.contacts
      results[[paste("class.nb.contacts",dist.class,sep=".")]]=cut(nb.contacts, breaks=c(1, 10, 20, 30, 40, max(nb.contacts)), include.lowest=T, labels=c("1-10", "11-20", "21-30", "30-40", ">40"))
      
      ## mean conservation score by gene
      
      mean.aln.score=tapply(filtered.contacts$align_score, factor(filtered.contacts$gene, levels=all.genes), mean, na.rm=T)
      
      results[[paste("mean.aln.score",dist.class,sep=".")]]=mean.aln.score
      results[[paste("class.aln.score", dist.class, sep=".")]]=cut(mean.aln.score, breaks=seq(from=0, to=1, length=6), include.lowest=T)

      ## synteny conservation

      filtered.synteny=synteny[which(synteny$id%in%filtered.contacts$id),]
      
      fr.cons.synt=tapply(filtered.synteny$target_dist,  factor(filtered.synteny$origin_gene, levels=all.genes), function(x) length(which(x<maxDistanceSyntenyTarget))/length(x))
      
      nb.enh.synt=as.numeric(table(factor(filtered.synteny$origin_gene, levels=all.genes)))
      
      ## if fewer than min.nb enhancers, we assign NA values to synteny conservation

      fr.cons.synt[which(nb.enh.synt < min.nb | nb.enh.synt > max.nb)]=NA

      results[[paste("fr.synteny.cons",dist.class,sep=".")]]=fr.cons.synt
      
      results[[paste("class.synteny.cons", dist.class, sep=".")]]=factor(fr.cons.synt<1,levels=c(TRUE, FALSE), labels=c("broken", "conserved"))
      
      ## contact conservation
      
      filtered.contact.cons=contact.cons[which(contact.cons$id%in%filtered.contacts$id),]

      fr.cons.contact=tapply(filtered.contact.cons$is_conserved, factor(filtered.contact.cons$origin_gene, levels=all.genes), function(x) length(which(x))/length(x))
      
      nb.enh.contact=as.numeric(table(factor(filtered.contact.cons$origin_gene, levels=all.genes)))
      
      ## if fewer than min.nb enhancers, we assign NA values to synteny conservation
      
      fr.cons.contact[which(nb.enh.contact < min.nb | nb.enh.contact > max.nb)]=NA
      
      results[[paste("fr.contact.cons",dist.class,sep=".")]]=fr.cons.contact
      
      results[[paste("class.contact.cons",dist.class, sep=".")]]=cut(fr.cons.contact, breaks=c(0, 0.2, 0.4, 0.6, 0.8,  1), include.lowest=T)
    }    
    
    results=as.data.frame(results)
    
    regland.conservation[[ref]][[enh]]=results
  }
}

#######################################################################################

save(regland.conservation, file=paste(pathFigures, "RData/data.regland.conservation.RData",sep=""))

#######################################################################################
#######################################################################################
