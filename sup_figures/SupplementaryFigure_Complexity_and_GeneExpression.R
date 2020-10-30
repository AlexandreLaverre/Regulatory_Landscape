######################################################################################################################
#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalData are defined based on the user name

sp="human"
sp_name="Human"

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/Fig6_", sp, "_common_cells.Rdata", sep=""))
load(paste(pathFigures, "RData/Fig6_", sp, "_all_samples_CardosoMoreira.Rdata", sep=""))

if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}
DivergenceMeasure = "CorrelationSpearman" # "EuclidianSimilarity" or "CorrelationSpearman" or SpearmanResidual

cells <- c("Bcell", "ESC", "adipo")
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

######################################################################################################################
xpos=seq(1, 5, 1)
smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]
xlim=c(0.5, 5.5)

CMPlot <- function(var, nb){
  if (var == "CorrelationSpearman"){ylim=c(0.55, 0.65)}else if (var == paste0("Tau", sp_name)){ylim=c(0.6, 0.75)}else{ylim=ylim=c(6, 10)}
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(class in levels(regland[[enh]]$class_nb_contact)){
      this.genes=rownames(regland[[enh]][which(regland[[enh]]$class_nb_contact == class),])
      
      xpos=seq(1, length(levels(regland[[enh]]$class_nb_contact)), 1)
      names(xpos) = levels(regland[[enh]]$class_nb_contact)
      x=xpos[class]+smallx[enh]
      
      b=boxplot(expdiv[[enh]][this.genes, var], plot=FALSE)
      med=median(expdiv[[enh]][this.genes, var])
      ci=as.numeric(b$conf)
      
      points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, ci[1], x, ci[2],  col=col.enhancers[enh])
    }
  }
  abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
  mtext(var, side=2, line=2.5, cex=0.9)
  
  if (nb == "3"){
    axis(side=1, cex.axis=1.2)
    mtext("Complexity", side=1, line=2, cex=0.9)
  }
}

CellTypesPlot <- function(var, nb){
  if (var == "Conservation"){ylim=c(0.2, 0.6)}else if (var == "MeanRPKM"){ylim=c(6, 35)}else{ylim=ylim=c(0, 1)}
  
  for (cell in cells){
    
    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
    
    if (var != ""){
      for(enh in enhancer.datasets[[sp]]){
        for(class in levels(data_cell[[cell]][[enh]]$class_nb_contact)){
          this.genes=rownames(data_cell[[cell]][[enh]][which(data_cell[[cell]][[enh]]$class_nb_contact == class),])
          
          xpos=seq(1, length(levels(data_cell[[cell]][[enh]]$class_nb_contact)), 1)
          names(xpos) = levels(data_cell[[cell]][[enh]]$class_nb_contact)
          x=xpos[class]+smallx[enh]
          
          b=boxplot(data_cell[[cell]][[enh]][this.genes, var], plot=FALSE)
          med=median(data_cell[[cell]][[enh]][this.genes, var])
          ci=as.numeric(b$conf)
          
          points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
          segments(x, ci[1], x, ci[2],  col=col.enhancers[enh])
        }
      }
      abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
      
      if (cell == "Bcell"){axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
        mtext(var, side=2, line=2.5, cex=0.9)}
      
      if (nb == "3"){axis(side=1, cex.axis=1.2)
        mtext("Complexity", side=1, line=2, cex=0.9)}
    }
    else{mtext(cell, side=1, line=1, cex=1.3)}
  }

}

MeasuresCM <- c(paste0("Tau", sp_name), "CorrelationSpearman", paste0(sp_name, "_MeanRPKM"))
MeasuresCellTypes <- c("", "Conservation", "MeanRPKM")

par(mfrow=c(3,4))
par(mai = c(0.4, 0.5, 0, 0.1)) # bottom, left, top, right

for (measure in 1:3){
  CMPlot(MeasuresCM[measure], measure)
  CellTypesPlot(MeasuresCellTypes[measure], measure)
}
