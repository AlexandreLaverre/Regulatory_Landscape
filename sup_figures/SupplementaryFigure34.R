######################################################################################################################
library(Hmisc)
#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  
  source("parameters.R")
}

##############################################################################
if(load){
  sp="human"
  load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
  load(paste(pathFigures, "RData/Fig6_", sp, "_SomaticOrgans_CardosoMoreira.Rdata", sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".gene.regland.conservation.RData", sep=""))
  if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}

  load=FALSE
}


################################################################################################################################
############################## Cardoso-Moreira  - Euclidean Similarity ##########################################################
expdiv$EuclideanSimilarity = 1-expdiv$EuclideanDistance

plot_profiles <- function(class_conserv, distances, xlab, xnames){
  smallx=c(-0.15, -0.075, 0.075, 0.15)
  names(smallx)=enhancer.datasets[[sp]]
  
  if (Measure == "corrected"){DivergenceMeasure = "CorrectedEuclideanSimilarity"; ylab="Residual Euclidean Similarity"; ylim=c(0, 0.02)
  }else{DivergenceMeasure = "EuclideanSimilarity"; ylab="Euclidean Similarity"; ylim=c(0.89, 0.94)}
  
  if (class_conserv == "class_cons_synt"){xmax=3}else{xmax=5}
  xlim=c(0.5, xmax+0.5)
  
  for (dist in distances){
    
    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
    
    for (enh in enhancer.datasets[[sp]]){
      regland = genes.conservation[[enh]][["obs"]][[dist]]
      genes = intersect(rownames(regland), rownames(expdiv))
      regland = regland[genes,]
      
      
      for (class in levels(regland[[class_conserv]])){
        this.genes=rownames(regland[which(regland[[class_conserv]] == class),])
        
        b=boxplot(expdiv[this.genes, DivergenceMeasure], plot=F)
        med=median(expdiv[this.genes, DivergenceMeasure])
        ci=as.numeric(b$conf)
        
        xpos=seq(1,  length(levels(regland[[class_conserv]])), 1)
        names(xpos) = levels(regland[[class_conserv]])
        
        x=xpos[class]+smallx[enh]
        
        points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
        segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
      }
    }
    
    abline(v=xpos[1:xmax-1]+0.5, lty=3, col="gray40")
    axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[class_conserv]]))), cex.axis=0.8)
    mtext(xnames, at=xpos, side=1, line=1, cex=0.8)
    mtext(xlab, side=1, line=2.5, cex=0.9)
    
    if (dist == "all"){
      axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
      mtext(ylab, side=2, line=2.5, cex=0.9)
      #mtext("d", side=3, at=0.45, font=2, cex=1.2, line=0.5)
      #legend("bottomleft", legend=enhancer.datasets[[sp]], pch=20,
      #       col=col.enhancers, cex=1, bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))
    }
    
    #mtext(dist, side=3, line=-1, cex=1)
  }
  return(expdiv)
}

################################################################################################################################
distances =  "all"  # c("25kb - 100kb", "100kb - 500kb", "500kb - 2Mb", "all")

pdf(file=paste(pathFigures, "SupplementaryFigure34_bis.pdf", sep=""), width = 8.5)

par(mfrow=c(2,3))
par(mai = c(0.5, 0.5, 0.5, 0.1)) # bottom, left, top, right

## Gene expression profiles uncorrected
Measure = "uncorrected"
expdiv <- plot_profiles("class_align_score", distances,  "Alignement score quantile", 1:5)
mtext("a", side=3, at=0.45, font=2, cex=1.2, line=0.5)
legend("bottomleft", legend=enhancer.datasets[[sp]], pch=20,
       col=col.enhancers, cex=1, bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))

expdiv <- plot_profiles("class_cons_synt", distances,  "Synteny Conservation", c("<75%", "75-99%", ">99%"))
expdiv <- plot_profiles("class_cons_cont", distances,  "Contact conservation", c("<1%", "25%", "50%", "75%", ">75%"))

## Gene expression profiles corrected
Measure = "corrected"
expdiv <- plot_profiles("class_align_score", distances,  "Alignement score quantile", 1:5)
mtext("b", side=3, at=0.45, font=2, cex=1.2, line=0.5)

expdiv <- plot_profiles("class_cons_synt", distances,  "Synteny Conservation", c("<75%", "75-99%", ">99%"))
expdiv <- plot_profiles("class_cons_cont", distances,  "Contact conservation", c("<1%", "25%", "50%", "75%", ">75%"))


dev.off()

