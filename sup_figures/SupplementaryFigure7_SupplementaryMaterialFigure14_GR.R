#########################################################################

source("../main_figures/parameters.R")

###########################################################################
  
## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

for(ref in c("human", "mouse")){
  ## load data
  load(paste(pathFigures, "RData/data.bootstrap.sequence.features.",ref,"RData",sep=""))
  
  tg=setdiff(c("human", "mouse"), ref)
  
  enh="ENCODE"

  if(ref=="human"){
    pdf(paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_S8.pdf", sep=""), width=6.85, height=6)
  } else{
    pdf(paste(pathFigures, "GenomeResearch_Figures/SupplementaryMaterialFigure14.pdf", sep=""), width=6.85, height=6)
  }
  
  m=matrix(rep(NA, 2*10), nrow=2)
  m[1,]=c(rep(1,2), rep(2, 4), rep(3, 4))
  m[2,]=c(rep(4,2), rep(5, 4), rep(6, 4))
  
  layout(m)
  
  ###########################################################################
  
  labels.genes=c("A", "D")
  labels.repeats=c("B", "E")
  labels.cons=c("C","F")
  
  names(labels.genes)=c("frag", "enh")
  names(labels.repeats)=c("frag", "enh")
  names(labels.cons)=c("frag", "enh")
  
  for(type in c("frag", "enh")){

    ## get objects
    stats.obs=get(paste(type, ".stats.obs", sep=""))
    stats.sim=get(paste(type, ".stats.sim", sep=""))
    
    mean.rep.obs=get(paste("mean.rep.",type,".obs",sep=""))
    ci.rep.low.obs=get(paste("ci.rep.",type,".low.obs",sep=""))
    ci.rep.high.obs=get(paste("ci.rep.",type,".high.obs",sep=""))
    
    mean.rep.sim=get(paste("mean.rep.",type,".sim",sep=""))
    ci.rep.low.sim=get(paste("ci.rep.",type,".low.sim",sep=""))
    ci.rep.high.sim=get(paste("ci.rep.",type,".high.sim",sep=""))
    
    mean.cons.obs=get(paste("mean.cons.repclass.",type,".obs", sep=""))
    ci.cons.low.obs=get(paste("ci.low.cons.repclass.",type,".obs", sep=""))
    ci.cons.high.obs=get(paste("ci.high.cons.repclass.",type,".obs", sep=""))
    
    mean.cons.sim=get(paste("mean.cons.repclass.",type,".sim", sep=""))
    ci.cons.low.sim=get(paste("ci.low.cons.repclass.",type,".sim", sep=""))
    ci.cons.high.sim=get(paste("ci.high.cons.repclass.",type,".sim", sep=""))
        
   ###########################################################################
    
    ## nb genes
    
    par(mar=c(6.75, 3.75, 2.1, 1.1))
    
    if(ref=="human"){
      ylim=range(c(0, 25))
    }
    
    if(ref=="mouse"){
      ylim=range(c(0, 40))
    } 
    
    plot(1, type="n", xlim=c(0.5, 2.5), ylim=ylim, xlab="", ylab="", axes=F)
    
    boxplot(stats.obs$nb_genes_500kb, col="white", border=dataset.colors[["Original"]], add=T, at=1, axes=F, outline=F, boxwex=0.5, notch=T)
    
    boxplot(stats.sim$nb_genes_500kb, pch=20, col="white", border=dataset.colors[["Simulated"]], add=T, at=2, axes=F, outline=F, boxwex=0.5, notch=T)
    
    axis(side=1, mgp=c(3, 0.5,0), at=1:2, labels=rep("", 2))
    axis(side=2, mgp=c(3, 0.75,0))
    
    mtext("number of genes within 500kb", side=2, line=2.5, cex=0.75)
    
    mtext(c("PCHi-C", "simulated"), at=1:2, side=1, line=0.75, cex=0.75, las=2)
    
    mtext(labels.genes[type], side=3, at=-0.75, font=2, line=1)
    
###########################################################################
    
    ## pc repeats
    
    par(mar=c(6.75, 3.5, 2.1, 1))
    
    xpos=1:length(mean.rep.obs)
    xlim=range(xpos)+c(-0.5, 0.5)
    
    ylim=range(c(ci.rep.low.obs, ci.rep.high.obs, ci.rep.low.sim, ci.rep.high.sim))
    smally=diff(ylim)/10
    ylim=ylim+c(-smally,smally)
    
    
    smallx=c(-0.1, 0.1)
    
    plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim)
    
    points(xpos+smallx[1], mean.rep.obs, col=dataset.colors["Original"], pch=20)
    segments(xpos+smallx[1], ci.rep.low.obs, xpos+smallx[1], ci.rep.high.obs, col=dataset.colors[["Original"]])
    
    points(xpos+smallx[2], mean.rep.sim, col=dataset.colors["Simulated"], pch=20)
    segments(xpos+smallx[2], ci.rep.low.sim, xpos+smallx[2], ci.rep.high.sim, col=dataset.colors[["Simulated"]])
    
    abline(v=xpos[-length(xpos)]+0.5, lty=3, col="gray40")
    
    axis(side=2, mgp=c(3, 0.75, 0), las=2)
    axis(side=1, mgp=c(3, 0.75, 0), at=xpos, labels=rep("", length(xpos)))
    
    mtext(names(mean.rep.obs), at=xpos, side=1, line=1, las=2, cex=0.75)
    
    mtext("number of genes within 500 kb", side=1, line=5, cex=0.75)
    mtext("% repetitive sequence", side=2, line=2.5, cex=0.75)
    
    mtext(labels.repeats[type], side=3, at=-1.15, font=2, line=1)
    
    if(type=="frag"){
      legend("topleft", col=dataset.colors, legend = c("PCHi-C data", "simulated data"), box.col="white", bg="white", pch=20, cex=1.1, inset=c(0.01,-0.15), xpd=NA)
    }
    
    ## conservation as a function of pc repeats
    
    par(mar=c(6.75, 3.5, 2.1, 1))
    
    xpos=1:length(mean.cons.obs)
    xlim=range(xpos)+c(-0.5, 0.5)
    
    ylim=range(c(ci.cons.low.obs, ci.cons.high.obs, ci.cons.low.sim, ci.cons.high.sim))
    smally=diff(ylim)/10
    ylim=ylim+c(-smally,smally)
    
    
    smallx=c(-0.1, 0.1)
    
    plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim)
    
    points(xpos+smallx[1], mean.cons.obs, col=dataset.colors["Original"], pch=20)
    segments(xpos+smallx[1], ci.cons.low.obs, xpos+smallx[1], ci.cons.high.obs, col=dataset.colors[["Original"]])
    
    points(xpos+smallx[2], mean.cons.sim, col=dataset.colors["Simulated"], pch=20)
    segments(xpos+smallx[2], ci.cons.low.sim, xpos+smallx[2], ci.cons.high.sim, col=dataset.colors[["Simulated"]])
    
    abline(v=xpos[-length(xpos)]+0.5, lty=3, col="gray40")
    
    axis(side=2, mgp=c(3, 0.75, 0), las=2)
    axis(side=1, mgp=c(3, 0.75, 0), at=xpos, labels=rep("", length(xpos)))
    
    mtext(names(mean.cons.obs), at=xpos, side=1, line=1, las=2, cex=0.75)
    
    mtext("repetitive sequence fraction", side=1, line=5, cex=0.75)
    mtext("% aligned sequence", side=2, line=2.5, cex=0.75)
    
    mtext(labels.cons[type], side=3, at=-0.7, font=2, line=1)
    
    if(type=="frag"){
      mtext("restriction fragments", side=4, line=-0.25, cex=0.75)
    } else{
      mtext("enhancers", side=4, line=-0.25, cex=0.75)
    }  
    
  }

  ## end of figure
  dev.off()
}
###########################################################################
