###########################################################################

path="/sps/biometr/necsulea/RegulatoryLandscapes/"
pathData=paste(path, "data/FANTOM5/", sep="")

###########################################################################

minpeaks=10000
mingenes=5000
minenhancers=10

minexpsamples.peaks=100
minexpsamples.enhancers=20
minexpsamples.genes=20

for(genome in c("mm9", "hg19")){
  sample.stats.peaks=read.table(paste(pathData, genome, "/SampleStatistics_MinTPM1_CAGE_Peaks.txt", sep=""), h=T, stringsAsFactors=F)
  sample.stats.enhancers=read.table(paste(pathData, genome, "/SampleStatistics_MinTPM1_Enhancers.txt", sep=""), h=T, stringsAsFactors=F)
  sample.stats.genes=read.table(paste(pathData, genome, "/SampleStatistics_MinTPM1_Genes.txt", sep=""), h=T, stringsAsFactors=F)
  
  peak.stats.peaks=read.table(paste(pathData, genome, "/PeakStatistics_MinTPM1_CAGE_Peaks.txt", sep=""), h=T, stringsAsFactors=F)
  peak.stats.enhancers=read.table(paste(pathData, genome, "/PeakStatistics_MinTPM1_Enhancers.txt", sep=""), h=T, stringsAsFactors=F)
  peak.stats.genes=read.table(paste(pathData, genome, "/PeakStatistics_MinTPM1_Genes.txt", sep=""), h=T, stringsAsFactors=F)

  par(mfrow=c(2,3))

  plot(density(sample.stats.peaks$NbPeaks), xlab="Nb expressed CAGE peaks per sample", ylab="Density", main="")
  plot(density(sample.stats.enhancers$NbPeaks), xlab="Nb expressed enhancers per sample", ylab="Density", main="")
  plot(density(sample.stats.genes$NbPeaks), xlab="Nb expressed genes per sample", ylab="Density", main="")
  
  plot(density(peak.stats.peaks$NbSamples), xlab="Nb samples with expression, CAGE peaks", ylab="Density", main="")
  plot(density(peak.stats.enhancers$NbSamples), xlab="Nb samples with expression, enhancers", ylab="Density", main="", xlim=c(0,200))
  plot(density(peak.stats.genes$NbSamples), xlab="Nb samples with expression, genes", ylab="Density", main="")
  
  ok.enhancers=peak.stats.enhancers$PeakID[which(peak.stats.enhancers$NbSamples>=minexpsamples.enhancers)]
  ok.peaks=peak.stats.peaks$PeakID[which(peak.stats.peaks$NbSamples>=minexpsamples.peaks)]
  ok.genes=peak.stats.genes$PeakID[which(peak.stats.genes$NbSamples>=minexpsamples.genes)]

  ok.samples1=sample.stats.peaks$SampleID[which(sample.stats.peaks$NbPeaks>=minpeaks)]
  ok.samples2=sample.stats.genes$SampleID[which(sample.stats.genes$NbPeaks>=mingenes)]
  ok.samples3=sample.stats.enhancers$SampleID[which(sample.stats.enhancers$NbPeaks>=minenhancers)]

  ok.samples=intersect(ok.samples3, intersect(ok.samples1, ok.samples2))
  
  writeLines(ok.enhancers, con=paste(pathData, genome, "/SelectedEnhancers.txt", sep=""))
  writeLines(ok.peaks, con=paste(pathData, genome, "/SelectedCAGEPeaks.txt", sep=""))
  writeLines(ok.genes, con=paste(pathData, genome, "/SelectedGenes.txt", sep=""))
  writeLines(ok.samples, con=paste(pathData, genome, "/SelectedSamples.txt", sep=""))
}

###########################################################################


