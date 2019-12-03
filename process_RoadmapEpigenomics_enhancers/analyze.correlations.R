###################################################################################

path="/beegfs/data/necsulea/RegulatoryLandscapes/"
pathResults=paste(path, "results/co_expression_analysis/RoadmapEpigenomics/", sep="")

###################################################################################

datasets=c("FOCS", "selected", "all")

for(dataset in datasets){

  print(dataset)

  ## real data
  real=read.table(paste(pathResults, dataset,"_expression_correlations_promoters_enhancers_in_contact_real_data.txt", sep=""), h=T, stringsAsFactors=F)

  real$PromoterStart=as.numeric(unlist(lapply(real$IDPromoter, function(x) unlist(strsplit(x, split=","))[2])))
  real$PromoterEnd=as.numeric(unlist(lapply(real$IDPromoter, function(x) unlist(strsplit(x, split=","))[3])))
  real$PromoterMidPos=(real$PromoterStart+real$PromoterEnd)/2
  
  real$EnhancerStart=as.numeric(unlist(lapply(real$IDEnhancer, function(x) unlist(strsplit(x, split=","))[2])))
  real$EnhancerEnd=as.numeric(unlist(lapply(real$IDEnhancer, function(x) unlist(strsplit(x, split=","))[3])))
  real$EnhancerMidPos=(real$EnhancerStart+real$EnhancerEnd)/2

  real$DistancePromoterEnhancer=abs(real$PromoterMidPos-real$EnhancerMidPos)

  real$DistanceClass=cut(real$DistancePromoterEnhancer, breaks=seq(from=25000, to=2.5e6, by=50000), include.lowest=T)

  ## simulated data
  
  simulated=read.table(paste(pathResults, dataset,"_expression_correlations_promoters_enhancers_in_contact_simulated_data.txt", sep=""), h=T, stringsAsFactors=F)

  simulated$PromoterStart=as.numeric(unlist(lapply(simulated$IDPromoter, function(x) unlist(strsplit(x, split=","))[2])))
  simulated$PromoterEnd=as.numeric(unlist(lapply(simulated$IDPromoter, function(x) unlist(strsplit(x, split=","))[3])))
  simulated$PromoterMidPos=(simulated$PromoterStart+simulated$PromoterEnd)/2
  
  simulated$EnhancerStart=as.numeric(unlist(lapply(simulated$IDEnhancer, function(x) unlist(strsplit(x, split=","))[2])))
  simulated$EnhancerEnd=as.numeric(unlist(lapply(simulated$IDEnhancer, function(x) unlist(strsplit(x, split=","))[3])))
  simulated$EnhancerMidPos=(simulated$EnhancerStart+simulated$EnhancerEnd)/2

  simulated$DistancePromoterEnhancer=abs(simulated$PromoterMidPos-simulated$EnhancerMidPos)

  simulated$DistanceClass=cut(simulated$DistancePromoterEnhancer, breaks=seq(from=25000, to=2.5e6, by=50000), include.lowest=T)


  pdf(file=paste("figures/PearsonCorrelation_",dataset,"_enhancers_promoters.pdf",sep=""), width=4.65, height=8.65)
  boxplot(real$PearsonCorrelation, simulated$PearsonCorrelation, notch=T, outline=F, names=c("observed", "simulated"), ylab="Pearson's correlation coefficient")
  dev.off()

  pdf(file=paste("figures/SpearmanCorrelation_",dataset,"_enhancers_promoters.pdf",sep=""), width=4.65, height=8.65)
  boxplot(real$SpearmanCorrelation, simulated$SpearmanCorrelation, notch=T, outline=F, names=c("observed", "simulated"), ylab="Spearman's correlation coefficient")
  dev.off()


  ## by distance
  
  medians.real=tapply(real$SpearmanCorrelation, real$DistanceClass, median, na.rm=T)
  ci.low.real=tapply(real$SpearmanCorrelation, real$DistanceClass, function(x) boxplot(x, plot=F, na.rm=T)$conf[1])
  ci.high.real=tapply(real$SpearmanCorrelation, real$DistanceClass, function(x) boxplot(x, plot=F, na.rm=T)$conf[2])

  medians.simulated=tapply(simulated$SpearmanCorrelation, simulated$DistanceClass, median, na.rm=T)
  ci.low.simulated=tapply(simulated$SpearmanCorrelation, simulated$DistanceClass, function(x) boxplot(x, plot=F, na.rm=T)$conf[1])
  ci.high.simulated=tapply(simulated$SpearmanCorrelation, simulated$DistanceClass, function(x) boxplot(x, plot=F, na.rm=T)$conf[2])

  ylim=range(c(ci.low.real, ci.high.real, ci.low.simulated, ci.low.simulated))

  pdf(file=paste("figures/SpearmanCorrelation_VariationWithDistance_",dataset,"_enhancers_promoters.pdf",sep=""), width=8.65, height=4.5)
  
  xpos=1:length(medians.real)
  plot(xpos, medians.real, type="b", col="red", xlim=c(0.5, max(xpos)+0.5), ylim=ylim, pch=20, xlab="", ylab="Spearman's correlation coefficient")
  mtext("Distance class (50 kb windows)", side=1, line=1, cex=0.8)
  segments(xpos, ci.low.real, xpos, ci.high.real, col="red")

  lines(xpos, medians.simulated, type="b", col="blue", pch=20)
  segments(xpos, ci.low.simulated, xpos, ci.high.simulated, col="blue")

  dev.off()
}

###################################################################################
