library(stringr)
library(data.table)

source("parameters.R")
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

pathCoverage="/home/laverre/Regulatory_landscape/result/HiCup/"

replicats = list()
coverage = list()

load(paste(pathFigures, "RData/data.fragment.statistics.RData",sep=""))

for(sp in c("mouse")){
  info=sampleinfo[[sp]]
  samples=info$Sample.ID
  author=word(info$Publication,1)
  first=T
  
  frag.obs <- fragment.statistics[[sp]][["original"]]
  frag.sim <- fragment.statistics[[sp]][["simulated"]]
  
  
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
  if (sp == "mouse"){maxlim=30000; breaks=5000; sumreads=50}else{maxlim=70000; breaks=5000; sumreads=100}
  hist(coverage[[sp]]$SumReads, breaks=breaks, xlim=c(0,maxlim), freq=F,
       xlab="Sum read counts by restriction fragment", ylab="Frequency", main="All fragments")
  
  low.covered.fragments=rownames(coverage[[sp]][which(coverage[[sp]]$SumReads < sumreads),])
  lower.20 = length(low.covered.fragments)
  #mtext(paste("Nb frag with sum reads <", sumreads, ":", lower.20, "on", nrow(coverage[[sp]])), line=-5)
  
  col = c(rgb(255,165,0,100, maxColorValue=255), rgb(0,0,128,100, maxColorValue=255))
  names(col) = c("obs", "sim")
  
  cov.obs = coverage[[sp]][frag.obs$ID,]$SumReads
  cov.sim = coverage[[sp]][frag.sim$ID,]$SumReads
  hist(cov.obs, col=col["obs"], xlim=c(0,30000), breaks=1000, freq=F,
       main="", xlab="Sum read counts by restriction fragment", ylab="Frequency", cex.lab=1.2)
  hist(cov.sim, col=col["sim"], breaks=1000, add=T, freq=F)
  legend("topright", legend=c("Observed restriction fragment", "Simulated"), fill=c(col["obs"], col["sim"]), bty='n')

  # cov.obs = coverage[[sp]][frag.obs$ID,]$MaxReads
  # cov.sim = coverage[[sp]][frag.sim$ID,]$MaxReads
  # hist(cov.obs, col=col["obs"], xlim=c(0,3000), breaks=3000, freq=F,
  #      main="", xlab="Max read counts by restriction fragment", ylab="Frequency", cex.lab=1.2)
  # hist(cov.sim, col=col["sim"], breaks=3000, add=T, freq=F)
  # legend("topright", legend=c("Observed restriction fragment", "Simulated"), fill=c(col["obs"], col["sim"]), bty='n')
  # 
  abline(v=147, col='red')
  dev.off()

  write.table(coverage[[sp]], file=paste(pathFinalData, "SupplementaryDataset1", sp, "reads.coverage.txt", sep="/"), quote=F, sep="\t")
  write.table(low.covered.fragments, file=paste(pathCoverage, sp, "low.covered.fragments.txt", sep="/"), quote=F, row.names = F, col.names = F)
  
  info$TotalReads = 0
  for (sample in info$Sample.ID){
    info[which(info$Sample.ID == sample),]$TotalReads = sum(coverage[[sp]][,grepl(sample, names(coverage[[sp]]))])
  }
  
  plot(info$TotalReads~info$Number.of.detected.interactions, xlab="Number of detected interactions", ylab="Total reads number")
  
}

