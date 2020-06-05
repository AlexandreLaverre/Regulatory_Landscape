############################################################################

pathExpression="../../../RegulatoryLandscapeManuscript/SupplementaryDataset6b/"

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

  exp=read.table(paste(pathExpression, sp, "/", filename, sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

}

############################################################################
