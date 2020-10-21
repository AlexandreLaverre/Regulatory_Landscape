############################################################################

## define paths

user=as.character(Sys.getenv()["USER"])

if(user=="laverre"){
  pathFinalData="/home/laverre/Manuscript/"
}

if(user=="necsulea"){
  pathFinalData="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
}


if(user=="ubuntu"){
  pathFinalData="/mnt/RegulatoryLandscapesManuscript/"
}

######################################################################################

pathExpression=paste(pathFinalData, "SupplementaryDataset6/", sep="")

############################################################################

compute.tau <- function(exp){
  if(max(exp)==0){
    return(NA)
  }
  
  n=length(exp)
  newexp=exp/max(exp)

  tau=sum(1-newexp)/(n-1)

  return(tau)
}

############################################################################

dataset="CardosoMoreira2019"
filename=paste("RPKMValues_",dataset,".txt", sep="")

############################################################################

for(sp in c("human", "mouse")){

  print(sp)

  exp=read.table(paste(pathExpression, sp, "/", filename, sep=""), h=T, stringsAsFactors=F, sep=" ", quote="\"")
  samples=colnames(exp)

  if(dataset=="CardosoMoreira2019"){
    tissue=unlist(lapply(samples, function(x) unlist(strsplit(x, split="\\."))[1]))
    age=unlist(lapply(samples, function(x) {y=unlist(strsplit(x, split="\\.")); return(paste(y[-c(1, length(y))], collapse="."))}))
    
    tissage=as.factor(paste(tissue, age, sep="."))
    
    avgexp=t(apply(exp, 1, function(x) tapply(as.numeric(x), tissage, mean)))
    
    avgexp=as.matrix(avgexp)

    write.table(avgexp, file=paste(pathExpression, sp, "/AverageRPKM_",dataset,".txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)
    
    tau.rpkm=apply(avgexp, 1, compute.tau)
    tau.logrpkm=apply(log2(avgexp+1), 1, compute.tau)

    overall.avg=apply(avgexp,1, mean)

    maxexp=apply(avgexp, 1, max)
    maxsample=apply(avgexp,1, function(x) colnames(avgexp)[which.max(x)])
    maxsample[which(maxexp==0)]=NA
    
    maxtiss=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="\\."))[1]))
    maxage=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="\\."))[2]))
    
    nbsamplesexp=apply(avgexp, 1, function(x) length(which(x>0)))
    
    results=data.frame("GeneID"=rownames(avgexp), "TauRPKM"=tau.rpkm, "TauLogRPKM"=tau.logrpkm, "MaxRPKM"=maxexp, "AverageRPKM"=overall.avg, "MaxSample"=maxsample, "MaxTissue"=maxtiss, "MaxAge"=maxage, "NbSamplesDetectedExpression"=nbsamplesexp, stringsAsFactors=F)
  }
  
  write.table(results, file=paste(pathExpression, sp, "/ExpressionStatistics_",dataset,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}

############################################################################
