library(stringr)
library(data.table)

source("parameters.R")
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

pathCoverage="/home/laverre/Regulatory_landscape/result/HiCup/"

replicats = list()
coverage = list()

for(sp in c("human")){
  info=sampleinfo[[sp]]
  samples=info$Sample.ID
  author=word(info$Publication,1)
  first=T
  
  for(i in 1:length(samples)){
    pathFiles=paste0(pathCoverage, "/", sp, "/", author[i], "/", samples[i], "/")
    files = list.files(pathFiles)
    
    for (replicat in files){
      replicat.name = str_remove(replicat, "_restriction_fragments_coverage.txt")
      print(replicat.name)
      replicats[[sp]]=c(replicats[[sp]], replicat.name)
      
      if (first){
        cov.sample = read.table(paste0(pathFiles, replicat), row.names = 4)
        coverage[[sp]] = cov.sample
        colnames(coverage[[sp]]) = c("chr", "start", "end", replicat.name)
        first=F
      
        }else{
        cov.sample = fread(paste0(pathFiles, replicat),  select=5)
        coverage[[sp]][[replicat.name]] = cov.sample$V5
        
        }
    }
  }
  
  coverage[[sp]]$SumReads=apply(coverage[[sp]][,replicats[[sp]]], 1, sum)
  coverage[[sp]]$MaxReads=apply(coverage[[sp]][,replicats[[sp]]], 1, max)
  
  pdf(paste(pathFigures, "CoverageDistribution.", sp, ".pdf", sep=""), width=6.85, height=7)
  par(mfrow=c(2,1))
  
  par(mai=c(1.2,0.8,0.3,0.5))
  boxplot(coverage[[sp]][,replicats[[sp]]], outline=F, notch=T, las=2, ylab="read counts", cex.axis=0.8)
  
  par(mai=c(0.8,0.8,0.3,0.5))
  if (sp == "mouse"){maxlim=20000; breaks=10000; sumreads=50}else{maxlim=70000; breaks=5000; sumreads=100}
  hist(coverage[[sp]]$SumReads, breaks=breaks, xlim=c(0,maxlim), xlab="Sum read counts by restriction fragment", ylab="Frequency", main="")
  
  low.covered.fragments=rownames(coverage[[sp]][which(coverage[[sp]]$SumReads < sumreads),])
  lower.20 = length(low.covered.fragments)
  mtext(paste("Nb frag with sum reads <", sumreads, ":", lower.20, "on", nrow(coverage[[sp]])), line=-5)
  dev.off()

  write.table(coverage[[sp]], file=paste(pathFinalData, "SupplementaryDataset1", sp, "reads.coverage.txt", sep="/"), quote=F, sep="\t")
  write.table(low.covered.fragments, file=paste(pathCoverage, sp, "low.covered.fragments.txt", sep="/"), quote=F, row.names = F, col.names = F)
  
  info$TotalReads = 0
  for (sample in info$Sample.ID){
    info[which(info$Sample.ID == sample),]$TotalReads = sum(coverage[[sp]][,grepl(sample, names(coverage[[sp]]))])
  }
  
  plot(info$TotalReads~info$Number.of.detected.interactions, xlab="Number of detected interactions", ylab="Total reads number")
  
}

