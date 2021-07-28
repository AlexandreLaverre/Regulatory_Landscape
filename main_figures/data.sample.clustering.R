#################################################################################

source("parameters.R")

#################################################################################

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))

minDistanceLongRange=250e3

subsample = FALSE
if (subsample){sample.suffix=".subsample"}else{sample.suffix=""}

#################################################################################

sample.clustering=list()

for(sp in c("human", "mouse")){
  info=sampleinfo[[sp]]
  samples=info$Sample.ID
  
  rownames(info)=samples
  
  observed=observed.contacts[[sp]]
  simulated=simulated.contacts[[sp]]
  
  # Sub-sampling to minimum number of interactions
  if (subsample){
    min.nb.interaction = min(apply(observed, 2, function(x) nrow(observed[which(!is.na(x)),])))
    
    for (sample in samples){
      treshold.value =  min(head(sort(observed[[sample]] ,decreasing=TRUE), n=min.nb.interaction))
      observed[which(observed[[sample]] < treshold.value), sample] = NA
    }
  }
  
  observed.longrange=observed[which(observed$distance>=minDistanceLongRange),]
  simulated.longrange=simulated[which(simulated$distance>=minDistanceLongRange),]
  
  ## matrix of similarity between pairs of samples, all distances

  mat.alldist.obs=matrix(rep(NA, length(samples)^2), nrow=length(samples))
  rownames(mat.alldist.obs)=samples
  colnames(mat.alldist.obs)=samples

  mat.alldist.sim=mat.alldist.obs

  ## only long-range distance
  mat.longrange.obs=mat.alldist.obs
  mat.longrange.sim=mat.alldist.obs
  
  for(i in 1:length(samples)){
    sample1=samples[i]
    
    mat.alldist.obs[sample1, sample1]=1
    mat.alldist.sim[sample1, sample1]=1
    mat.longrange.obs[sample1, sample1]=1
    mat.longrange.sim[sample1, sample1]=1
    
    if(i<length(samples)){
      for(j in (i+1):length(samples)){
        sample2=samples[j]

        print(paste(i, j, sample1, sample2))

        mat.alldist.obs[sample1, sample2]=length(which((!is.na(observed[,sample1])) & (!is.na(observed[,sample2]))))/length(which((!is.na(observed[,sample1])) | (!is.na(observed[,sample2]))))
        mat.alldist.sim[sample1, sample2]=length(which((!is.na(simulated[,sample1])) & (!is.na(simulated[,sample2]))))/length(which((!is.na(simulated[,sample1])) | (!is.na(simulated[,sample2]))))

        mat.longrange.obs[sample1, sample2]=length(which((!is.na(observed.longrange[,sample1])) & (!is.na(observed.longrange[,sample2]))))/length(which((!is.na(observed.longrange[,sample1])) | (!is.na(observed.longrange[,sample2]))))
        mat.longrange.sim[sample1, sample2]=length(which((!is.na(simulated.longrange[,sample1])) & (!is.na(simulated.longrange[,sample2]))))/length(which((!is.na(simulated.longrange[,sample1])) | (!is.na(simulated.longrange[,sample2]))))


        mat.alldist.obs[sample2, sample1]=mat.alldist.obs[sample1, sample2]
        mat.alldist.sim[sample2, sample1]=mat.alldist.sim[sample1, sample2]
        mat.longrange.obs[sample2, sample1]=mat.longrange.obs[sample1, sample2]
        mat.longrange.sim[sample2, sample1]=mat.longrange.sim[sample1, sample2]
                         
      }
    }
  }

  diff.alldist=mat.alldist.obs-mat.alldist.sim
  diff.longrange=mat.longrange.obs-mat.longrange.sim

  hclust.obs=hclust(as.dist(1-mat.alldist.obs))
  hclust.alldist=hclust(as.dist(1-diff.alldist))
  hclust.longrange=hclust(as.dist(1-diff.longrange))

  sample.order.obs=rownames(mat.alldist.obs)[order.dendrogram(as.dendrogram(hclust.obs))]
  sample.order.alldist=rownames(diff.alldist)[order.dendrogram(as.dendrogram(hclust.alldist))]
  sample.order.longrange=rownames(diff.longrange)[order.dendrogram(as.dendrogram(hclust.longrange))]

  ## store results

  sample.clustering[[sp]]=list("mat.alldist.obs"=mat.alldist.obs, "mat.alldist.sim"=mat.alldist.sim,
                               "mat.longrange.obs"=mat.longrange.obs, "mat.longrange.sim"=mat.longrange.sim,
                               "hclust.obs"=hclust.obs, "hclust.alldist"=hclust.alldist, "hclust.longrange"=hclust.longrange, 
                               "sample.order.obs"=sample.order.obs, "sample.order.alldist"=sample.order.alldist, "sample.order.longrange"=sample.order.longrange)
}

#################################################################################

save(sample.clustering, file=paste(pathFigures, "RData/data.sample.clustering", sample.suffix, ".RData", sep=""))

#################################################################################
