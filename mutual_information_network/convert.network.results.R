###########################################################################

path="/sps/biometr/necsulea/RegulatoryLandscapes/"
pathMI=paste(path, "results/mutual_information_network/", sep="")
pathFANTOM=paste(path, "data/FANTOM5/", sep="")
pathAnnot=paste(path, "data/ensembl_annotations/", sep="")

###########################################################################

ensrelease=94

###########################################################################

species=c("Mouse", "Human")
names(species)=c("mm9", "hg19")

newassemblies=c("mm10", "hg38")
names(newassemblies)=c("mm9", "hg19")

###########################################################################

for(genome in c("hg19")){
  newassembly=newassemblies[genome]
  sp=species[genome]

  ## enhancer coordinates lifted to new assembly

  enhancer.coords=read.table(paste(pathFANTOM, genome,"/", tolower(sp), "_lifted_to_", newassembly,"_permissive_enhancers_phase_1_and_2.bed", sep=""), h=F, stringsAsFactors=F, sep="\t")

  rownames(enhancer.coords)=enhancer.coords$V4
  colnames(enhancer.coords)[1]="chr"
  colnames(enhancer.coords)[2]="start"
  colnames(enhancer.coords)[3]="end"

  ## gene names
  gene.info=read.table(paste(pathAnnot, sp, "/GeneInfo_Ensembl",ensrelease, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  weirdchr=unique(gene.info$name[c(grep("PATCH", gene.info$name), grep("HSCHR", gene.info$name))])
  gene.info=gene.info[which(!gene.info$name%in%weirdchr),]
  rownames(gene.info)=gene.info$stable_id
  
  ## gene names

  gene.names=read.table(paste(pathAnnot, sp, "/GeneNames_Ensembl",ensrelease, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  colnames(gene.names)=c("GeneID", "GeneName")
  gene.names=gene.names[which(gene.names$GeneID%in%gene.info$stable_id),] ## remove genes on weird chromosomes first - they may share names with genes on regular chromosomes
  
  dupli.names=gene.names$GeneName[which(duplicated(gene.names$GeneName))] ## remove ambiguous names
  gene.names=gene.names[which(!gene.names$GeneName%in%dupli.names),]
  rownames(gene.names)=gene.names$GeneName
 
  ## mutual information network
  
  for(method in c("withBonferroni", "withoutBonferroni")){
    print(paste("reading network for ", method, sep=""))
    
    network=read.table(paste(pathMI, genome, "/aracne_replicates/network_", method, ".txt", sep=""), h=T,stringsAsFactors=F)
    colnames(network)[1]="EnhancerID"
    colnames(network)[2]="GeneName"

    ## we can also have enhancer-enhancer interactions, we don't keep them here
    ## select only enhancers that were projected on new gene assembly and genes with unambiguous names
    
    network=network[which(network$EnhancerID%in%rownames(enhancer.coords) & network$GeneName%in%rownames(gene.names)),]
    
    network$EnhancerChr=enhancer.coords[network$EnhancerID, "chr"]
    network$EnhancerStart=enhancer.coords[network$EnhancerID, "start"]
    network$EnhancerEnd=enhancer.coords[network$EnhancerID, "end"]
    
    network$EnhancerChr=unlist(lapply(network$EnhancerChr, function(x) substr(x, 4, nchar(x))))
    
    network$GeneID=gene.names[network$GeneName, "GeneID"]
    network$GeneChr=gene.info[network$GeneID, "name"]
    network$GeneStart=gene.info[network$GeneID, "seq_region_start"]
    network$GeneEnd=gene.info[network$GeneID, "seq_region_end"]
    network$GeneStrand=gene.info[network$GeneID, "seq_region_strand"]
    network$GeneType=gene.info[network$GeneID, "biotype"]
    
    network$InteractionType=rep("trans", dim(network)[1])
    network$InteractionType[which(network$GeneChr==network$EnhancerChr)]="cis"

    maxstart=apply(network[,c("EnhancerStart", "GeneStart")], 1, max)
    minend=apply(network[,c("EnhancerEnd", "GeneEnd")], 1, min)

    network$OverlapCoordinates=rep("no", dim(network)[1])
    network$OverlapCoordinates[which(network$InteractionType=="cis" & maxstart<=minend)]="yes"

    write.table(network, file=paste(pathMI, genome, "/NetworkInfo_",method,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
   
  }
}

###########################################################################

