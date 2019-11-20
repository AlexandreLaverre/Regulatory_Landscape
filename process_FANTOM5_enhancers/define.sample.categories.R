###########################################################################

path="/sps/biometr/necsulea/RegulatoryLandscapes/"
pathData=paste(path, "data/FANTOM5/", sep="")

###########################################################################

forbidden=c("whole body", "unclassified", "cancer", "leukemia", "carcinoid", "tumor", "cell line",  "treated", "induction", "response", "nuclear fraction", "ribopure", "pluriselect", "infection", "Universal")

###########################################################################

for(genome in c("mm9", "hg19")){

  ## statistics for enhancers and genes

  stats.enhancers=read.table(paste(pathData, genome, "/SampleStatistics_MinTPM1_Enhancers.txt", sep=""), h=T, stringsAsFactors=F)
  stats.genes=read.table(paste(pathData, genome, "/SampleStatistics_MinTPM1_Genes.txt", sep=""), h=T, stringsAsFactors=F)

  common.samples=intersect(stats.enhancers$SampleID, stats.genes$SampleID)
  
  sampleinfo=readLines(paste(pathData, genome, "/SampleDescription.txt", sep=""))
  sampleid=unlist(lapply(sampleinfo, function(x) grep("^CN", unlist(strsplit(x, split="\\.")), value=T)))
  description=unlist(lapply(sampleinfo, function(x) unlist(strsplit(x, split="\\."))[1]))
  description=unlist(lapply(description, function(x) unlist(strsplit(x, split="tpm of "))[2]))
  
  results=data.frame("SampleID"=sampleid, "Description"=description, stringsAsFactors=F)
  results=results[which(results$SampleID%in%common.samples),]

  notok=unique(results$SampleID[unlist(lapply(forbidden, function(x) grep(x, results$Description)))])
  print(paste("removing", length(notok), "samples with forbidden keywords"))

  results=results[which(!results$SampleID%in%notok),]

  results$ReplicateInfo=rep(NA, dim(results)[1])
  results$ReplicateInfo[grep(", biol_rep", results$Description)]=", biol_rep"
  results$ReplicateInfo[grep(", biol_ rep", results$Description)]=", biol_ rep"
  results$ReplicateInfo[grep(", pool", results$Description)]=", pool"
  results$ReplicateInfo[grep(", donor", results$Description)]=", donor"
  results$ReplicateInfo[grep(" donor", results$Description)]=" donor"
  results$ReplicateInfo[grep(", rep", results$Description)]=", rep"
  results$ReplicateInfo[grep(":", results$Description)]=":"

  print(paste("NA replicate info", length(which(is.na(results$ReplicateInfo)))))
  print(table(results$ReplicateInfo))

  results$SampleType=rep(NA, dim(results)[1])
  w=which(!is.na(results$ReplicateInfo))
  results$SampleType[w]=unlist(lapply(w, function(x) unlist(strsplit(results$Description[x], split=results$ReplicateInfo[x]))[1]))

  results$SampleType[which(is.na(results$ReplicateInfo))]=results$Description[which(is.na(results$ReplicateInfo))]

  print(paste(length(unique(results$SampleType)), "unique sample types"))
  print(paste(length(which(is.na(results$SampleType))), "undefined sample types"))

  results$TypeID=paste("ST",as.numeric(as.factor(results$SampleType)), sep="")
  
  write.table(results, file=paste(pathData, genome, "/SampleTypes.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

  
}

###########################################################################
