#######################################################################################

library(data.table)

options(stringsAsFactors = FALSE)

source("parameters.R")

#######################################################################################

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData", sep=""))

## we load sequence and synteny conservation later

#######################################################################################

## minDistance and maxDistance defined in parameters.R

dist.classes=c("all", "shortrange", "longrange")
min.distances=c(minDistance, minDistance, 100001)
max.distances=c(maxDistance, 100000, maxDistance)

## min number of enhancers for contacts and synteny conservation

min.nb=5

#######################################################################################

regland.conservation=list()

for(ref in c("human", "mouse")){
  
  tg=setdiff(c("human", "mouse"), ref)

  sampleinfo.ref=sampleinfo[[ref]]
  sampleinfo.tg=sampleinfo[[tg]]

  samples.ref=sampleinfo.ref[,"Sample.ID"]
  samples.tg=sampleinfo.tg[,"Sample.ID"]

  regland.conservation[[ref]]=list()
    
  for(enh in enhancer.datasets[[ref]]){
    print(paste(ref, enh))

    ## we start from the gene-enhancer contact data

    contacts=gene.enhancer.contacts[[ref]][[enh]][["real"]]
 
    ## we select only previously filtered orthologous genes

    contacts=contacts[which(contacts$gene%in%ortho[,ref]),]

    ## we add sequence conservation column

    load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref,"2", tg,".RData", sep=""))
    contacts$align_score=pcungapped[contacts$enhancer]

    ## we select enhancers that have an actual sequence conservation score, no NA values

    contacts=contacts[which(!is.na(contacts$align_score)),]

    ## we load synteny conservation - already filtered, min sequence conservation >= 10% threshold
    
    load(paste(pathFigures, "/RData/data.synteny.conservation.", ref, ".RData", sep=""))

    synteny=conserv_synteny[[enh]][[tg]][["synt_obs"]]

    ## we load contact conservation - already filtered, min sequence conservation 0.4

    contact.cons=contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]][["obs"]]

    contact.cons$is_conserved=apply(contact.cons[,samples.tg], 1, function(x) any(x>0))
    
    ## full gene list: ortho genes in contact dataset

    all.genes=unique(contacts$gene)

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

      results[[paste("nb.contacts",dist.class,sep=".")]]=nb.contacts

      ## median conservation score by gene

      median.aln.score=tapply(filtered.contacts$align_score, factor(filtered.contacts$gene, levels=all.genes), median, na.rm=T)
      
      results[[paste("median.aln.score",dist.class,sep=".")]]=median.aln.score

      ## synteny conservation

      filtered.synteny=synteny[which(synteny$id%in%filtered.contacts$id),]

      fr.cons.synt=tapply(filtered.synteny$target_dist,  factor(filtered.synteny$origin_gene, levels=all.genes), function(x) length(which(x<maxDistanceSyntenyTarget))/length(x))

      nb.enh.synt=as.numeric(table(factor(filtered.synteny$origin_gene, levels=all.genes)))

      ## if fewer than min.nb enhancers, we assign NA values to synteny conservation

      fr.cons.synt[which(nb.enh.synt<min.nb)]=NA

      results[[paste("fr.synteny.cons",dist.class,sep=".")]]=fr.cons.synt

      ## contact conservation

      filtered.contact.cons=contact.cons[which(contact.cons$id%in%filtered.contacts$id),]

      fr.cons.contact=tapply(filtered.contact.cons$is_conserved, factor(filtered.contact.cons$origin_gene, levels=all.genes), function(x) length(which(x))/length(x))
      
      nb.enh.contact=as.numeric(table(factor(filtered.contact.cons$origin_gene, levels=all.genes)))
      
      ## if fewer than min.nb enhancers, we assign NA values to synteny conservation

      fr.cons.contact[which(nb.enh.contact<min.nb)]=NA
      
      results[[paste("fr.contact.cons",dist.class,sep=".")]]=fr.cons.contact
    }

    results=as.data.frame(results)

    regland.conservation[[ref]]=results
  }
}

#######################################################################################

save(regland.conservation, file=paste(pathFigures, "RData/data.", ref, ".regland.conservation.RData",sep=""))

#######################################################################################
#######################################################################################
