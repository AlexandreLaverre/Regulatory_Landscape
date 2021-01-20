############################################################################

path="/beegfs/data/necsulea/RegulatoryLandscapes/"
pathFOCS=paste(path, "data/FOCS/", sep="")
pathResults=paste(path, "results/co_expression_analysis/", sep="")

options(stringsAsFactors=F)
options(digits=2) ## to make files lighter

########################################################################

datasets=list()
datasets[["human"]]=c("ENCODE", "FANTOM5", "GRO-seq", "RoadmapEpigenomics")
datasets[["mouse"]]=c("FANTOM5")

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


for(sp in c("human", "mouse")){
  for(dataset in datasets[[sp]]){
    print(paste(sp, dataset))
    
    exp=read.table(paste(pathFOCS, sp, "/", dataset, "/enhancer_coords_activity.txt", sep=""),  h=T, stringsAsFactors=F)
 
    samples=colnames(exp)[5:dim(exp)[2]] ## first 4 columns are coordinates
        
    tau.tpm=apply(exp[,samples], 1, compute.tau)
    tau.logtpm=apply(log2(exp[,samples]+1), 1, compute.tau)

    overall.avg=apply(exp[,samples],1, mean)

    maxexp=apply(exp[,samples], 1, max)
   
    nbsamplesexp=apply(exp[,samples], 1, function(x) length(which(x>0)))

     results=data.frame("id"=exp$id, "chr"=exp$chr, "start"=exp$start, "end"=exp$end, "Tau"=tau.tpm, "TauLog"=tau.logtpm, "MaxExp"=maxexp, "AverageExp"=overall.avg, "NbSamplesExp"=nbsamplesexp,stringsAsFactors=F)

    for(minTPM in c(0.1, 0.5, 1)){
      nbsampleshighexp=apply(exp[,samples], 1, function(x) length(which(x>=minTPM)))
      
      results[,paste("NbSamplesTPM",minTPM, sep="")]=nbsampleshighexp
    }
    
    write.table(results, file=paste(pathResults, sp, "/", dataset, "/enhancer_activity_statistics.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
    
  }
}

############################################################################
