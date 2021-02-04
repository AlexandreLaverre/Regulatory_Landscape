###########################################################################

library(data.table)

source("parameters.R")

pathContacts=paste(pathFinalData, "SupplementaryDataset4/", sep="")

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))

###########################################################################

gene.enhancer.contacts=list()

for(sp in c("human", "mouse")){
  print(sp)

  ## gene annotations 

  annot=gene.annot[[sp]]
  pc.genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
  
  ## previously filtered fragment contacts

  frag.contact.real=observed.contacts[[sp]]
  frag.contact.sim=simulated.contacts[[sp]]

  frag.contact.real$idcontact=paste(frag.contact.real$id_bait, frag.contact.real$id_frag, sep="-")
  frag.contact.sim$idcontact=paste(frag.contact.sim$id_bait, frag.contact.sim$id_frag, sep="-")

  ## sample info

  info=sampleinfo[[sp]]
  rownames(info)=info$Sample.ID
  
  samples=info$Sample.ID 
  celltypes=info$Broad.cell.type.or.tissue
  names(celltypes)=samples

  ## gene enhancer contacts
  
  gene.enhancer.contacts[[sp]]=list()
  
  for(enh in enhancer.datasets[[sp]]){
    print(enh)
    
    ## previously filtered enhancers
    
    enhancers.real=enhancer.statistics[[sp]][[enh]][["original"]]
    enhancers.sim=enhancer.statistics[[sp]][[enh]][["simulated"]]

    ## read contact data
    
    real <- fread(paste(pathContacts, sp, "/", enh, "/gene_enhancer_contacts_original_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    class(real) <- "data.frame"

    real$id=paste(real$gene, real$enhancer, sep="-")
        
    print(paste(nrow(real), "contacts before filtering, observed data"))
    
    sim<-fread(paste(pathContacts, sp, "/", enh, "/gene_enhancer_contacts_simulated_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    class(sim) <- "data.frame"

    sim$id=paste(sim$gene, sim$enhancer, sep="-")
    
    print(paste(nrow(sim), "contacts before filtering, simulated data"))

    ## take only previously filtered enhancers

    real=real[which(real$enhancer%in%enhancers.real$enh),]
    sim=sim[which(sim$enhancer%in%enhancers.sim$enh),]

    ## expand tables and select interactions for which the bait and the contacted fragments are in the previously filtered contact list
    ## real data

    list.baits.real=lapply(real$baits, function(x) unlist(strsplit(x, split=",")))
    list.frags.real=lapply(real$contacts, function(x) unlist(strsplit(x, split=",")))

    real$nb_baits=unlist(lapply(list.baits.real, length))
    real$nb_frag=unlist(lapply(list.frags.real, length))

    real$nb_pairs=real$nb_baits*real$nb_frag
    
    real.gene=rep(real$gene, real$nb_pairs)
    real.enhancer=rep(real$enhancer, real$nb_pairs)
   
    real.baits=unlist(lapply(1:nrow(real), function(x) rep(list.baits.real[[x]], length(list.frags.real[[x]]))))
    real.frags=unlist(lapply(1:nrow(real), function(x) rep(list.frags.real[[x]], each=length(list.baits.real[[x]]))))

    real.expanded=data.frame("gene"=real.gene,"enhancer"=real.enhancer, "bait"=real.baits, "frag"=real.frags, stringsAsFactors=F)
    real.expanded$idcontact=paste(real.expanded$bait, real.expanded$frag, sep="-")

    real.expanded=real.expanded[which(real.expanded$idcontact%in%frag.contact.real$idcontact),]

    real.expanded$id=paste(real.expanded$gene, real.expanded$enhancer, sep="-")
    real=real[which(real$id%in%real.expanded$id),]

    print(paste(nrow(real), "contacts after filtering, observed data"))
    
    ## simulated data
    list.baits.sim=lapply(sim$baits, function(x) unlist(strsplit(x, split=",")))
    list.frags.sim=lapply(sim$contacts, function(x) unlist(strsplit(x, split=",")))
    
    sim$nb_baits=unlist(lapply(list.baits.sim, length))
    sim$nb_frag=unlist(lapply(list.frags.sim, length))
    
    sim$nb_pairs=sim$nb_baits*sim$nb_frag
    
    sim.gene=rep(sim$gene, sim$nb_pairs)
    sim.enhancer=rep(sim$enhancer, sim$nb_pairs)
   
    sim.baits=unlist(lapply(1:nrow(sim), function(x) rep(list.baits.sim[[x]], length(list.frags.sim[[x]]))))
    sim.frags=unlist(lapply(1:nrow(sim), function(x) rep(list.frags.sim[[x]], each=length(list.baits.sim[[x]]))))
  
    sim.expanded=data.frame("gene"=sim.gene,"enhancer"=sim.enhancer, "bait"=sim.baits, "frag"=sim.frags, stringsAsFactors=F)
    sim.expanded$idcontact=paste(sim.expanded$bait, sim.expanded$frag, sep="-")

    sim.expanded=sim.expanded[which(sim.expanded$idcontact%in%frag.contact.sim$idcontact),]

    sim.expanded$id=paste(sim.expanded$gene, sim.expanded$enhancer, sep="-")
    sim=sim[which(sim$id%in%sim.expanded$id),]

    print(paste(nrow(sim), "contacts after filtering, simulated data"))

    ## number of cells in which there are interactions
    
    real$nb_cell <- apply(real[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
    sim$nb_cell <- apply(sim[,samples], 1, function(x) length(unique(celltypes[which(!is.na(x))])))
    
    ## select interactions in the accepted distance range
    real=real[which(real$dist>=minDistance & real$dist<=maxDistance),]
    sim=sim[which(sim$dist>=minDistance & sim$dist<=maxDistance),]

    ## select protein-coding genes

    real=real[which(real$gene%in%pc.genes),]
    sim=sim[which(sim$gene%in%pc.genes),]

    ## save data
    
    gene.enhancer.contacts[[sp]][[enh]]=list("real"=real, "simulated"=sim)
  }
}

###########################################################################

save(list=c("gene.enhancer.contacts"), file=paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))

###########################################################################


