#############################################################################################

source("../main_figures/parameters.R")

load(paste(pathFigures, "RData/data.bait.euclidean.distances.RData",sep=""))

######################################################################################
pdf(file=paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_X_euclidean.distance.by.bait.pdf", sep=""), width=6.85, height=7)

m=matrix(rep(NA, 3*3), nrow=3)

m[1,]=1:3
m[2,]=4:6
m[3,]=7:9

layout(m)

######################################################################################
label=1

for(sp in c("human", "mouse")){
  if (sp == "human"){sample="Mac0"}else{sample="preB_aged"}
  
  ######################################################################################
  ## Plot some example
  par(mar=c(3.1, 3.1, 2.1, 1.1))
  
  plot(example[[sp]]$distrib.sim, type="h", xlab="", ylab="", main="example bait", axes=F, col=dataset.colors["Original"])
  lines(example[[sp]]$distrib.obs, col=dataset.colors["Simulated"])
  lines(distrib.all[[sp]]$density, col="red", cex=1.2)
  
  legend("topright", col=c(dataset.colors,"red"),  bty="n", lty=1,
         legend=c(paste0("PCHi-C, N=",example[[sp]]$Nobs),paste0("simulated, N=",example[[sp]]$Nsim), "average empirical,"))
  
  mtext(paste0("N=",distrib.all[[sp]]$N), side=3, line=-4.5, at=280, cex=0.7)
  
  at = c("human"= c(0.05, 0.04), "mouse"=c(0.038, 0.03))
  mtext(paste0("eucl. dist.=", signif(example[[sp]]$dist.obs, digits=2)), side=4, cex=0.75, las=2, line=-7, at=at[paste0(sp, 1)], col=dataset.colors["Original"]) 
  mtext(paste0("eucl. dist.=", signif(example[[sp]]$dist.sim, digits=2)), side=4, cex=0.75, las=2, line=-7, at=at[paste0(sp, 2)], col=dataset.colors["Simulated"]) 
  
  xaxs=c(0, 100, 200, 300, 400)
  xaxslabels=c("0", "0.5", "1", "1.5", "2")  
  axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.95)
  mtext("distance bait-fragment (Mb)", side=1, line=1.5, cex=0.75)
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.95)
  mtext("density of nb contacts", side=2, line=2, cex=0.75)
  
  mtext(LETTERS[label], side=3, at=-50, cex=1, font=2, line=0.8)
  label = label+1
  
  ######################################################################################
  # Bait obs vs all obs
  med=signif(mean(distance.bait.obs.all[[sp]][[sample]]), digits=2)
  sd=signif(sd(distance.bait.obs.all[[sp]][[sample]]), digits=2)
  
  hist(distance.bait.obs.all[[sp]][[sample]], breaks=25, xlab="", ylab="", main="", cex.lab=1.2, xlim=c(0,0.8), col=dataset.colors["Original"])

  mtext(paste("mean = ", med, "\n"," sd = ", sd, sep=""), side=3, line=-1.5, at=0.6, cex=0.75)

  mtext("euclidean dist. to average empirical", side=1, line=1.85, cex=0.75)
  mtext("nb. baits", side=2, line=2.1, cex=0.75)
  
  mtext(LETTERS[label], side=3, at=-0.1, cex=1.1, font=2, line=0.85)
  label = label+1
  
  # Bait sim vs all sim
  med=signif(mean(distance.bait.sim.all[[sp]][[sample]]), digits=2)
  sd=signif(sd(distance.bait.sim.all[[sp]][[sample]]), digits=2)
  
  hist(distance.bait.sim.all[[sp]][[sample]], breaks=25, xlab="", ylab="", main="", cex.lab=1.2, xlim=c(0,0.8), col=dataset.colors["Simulated"])
  
  mtext(paste("mean = ", med, "\n"," sd = ", sd, sep=""), side=3, line=-1.5, at=0.6, cex=0.75)
  mtext(paste(sp,"\n",sample, "\n", "N=", length(distance.bait.obs.all[[sp]][[sample]]),  sep=""), side=4, cex=0.75, las=2, line=-3.7) 
  mtext("euclidean dist. to average empirical", side=1, line=1.85, cex=0.75)
  mtext("nb. baits", side=2, line=2.1, cex=0.75)
  
  mtext(LETTERS[label], side=3, at=-0.1, cex=1.1, font=2, line=0.85)
  label = label+1
  
  ######################################################################################
  
}

for(sp in c("human", "mouse")){
  ## mean distance, human
  par(mar=c(3.1, 3.1, 2.1, 1.1))
  
  plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=c(0.17, 0.42), ylim=c(0.17, 0.42))
  
  # confidence interval for PCHIC
  segments(x0=mean.dist.bait.obs[[sp]]-sd.dist.bait.obs[[sp]], col=dataset.colors["Original"],
           x1=mean.dist.bait.obs[[sp]]+sd.dist.bait.obs[[sp]], y0=mean.dist.bait.sim[[sp]], lwd=0.05)
  
  # confidence interval for simul
  segments(y0=mean.dist.bait.sim[[sp]]-sd.dist.bait.sim[[sp]], col=dataset.colors["Simulated"],
           y1=mean.dist.bait.sim[[sp]]+sd.dist.bait.sim[[sp]], x0=mean.dist.bait.obs[[sp]], lwd=0.05)
  
  points(mean.dist.bait.obs[[sp]], mean.dist.bait.sim[[sp]], pch=20, cex=0.6)

  abline(0,1)
  
  axis(side=1, mgp=c(3, 0.5, 0))
  axis(side=2, mgp=c(3, 0.5, 0))
  
  box()
  
  mtext("mean euclidean dist., PCHi-C data", side=1, line=1.75, cex=0.75)
  mtext("mean euclidean dist., simulated", side=2, line=2.1, cex=0.75)
  
  mtext(LETTERS[label], side=3, at=0.14, cex=1.1, font=2, line=0.85)
  label = label+1
  
}

dev.off()

# Bait obs vs bait simul

# N.baits=length(distance.bait.obs.sim)
# med=signif(mean(distance.bait.obs.sim), digits=2)
# sd=signif(sd(distance.bait.obs.sim), digits=2)
# 
# hist(distance.bait.obs.sim, breaks=100, xlab="Euclidean Distance", main="Bait-obs vs Bait-sim", cex.lab=1.2)
# legend("right", legend=c(paste0("mean = ", med), paste0("sd = ", sd)), bty='n')

