######################################################################################################################
setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalData are defined based on the user name

sp="human"
sp_name="Human"

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/Fig6_", sp, "_common_cells.Rdata", sep=""))
load(paste(pathFigures, "RData/Fig6_", sp, "_SomaticOrgans_CardosoMoreira.Rdata", sep=""))
expdiv$EuclideanSimilarity <- 1-expdiv$EuclideanDistance

if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}

cells <- c("Bcell", "ESC", "adipo")
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

#############################################################################################################
###################################### Expression profiles ##################################################
xpos=seq(1, 5, 1)
smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]
xlim=c(0.5, 5.5)

CMPlot <- function(var, plot.nb){
  if (var == "ResidualEuclideanSimilarity"){ylim=c(0.01, 0.04)
  }else if (var == "ResidualSpearman"){ylim=c(0.02, 0.1)
  }else if (var == "CorrectedSpearman"){ylim=c(0.02, 0.1)
  }else if (var == "CorrectedEuclideanSimilarity"){ylim=c(0, 0.02)
  }else if (var == "EuclideanSimilarity"){ylim=c(0.9, 0.95)
  }else if (var == "CorrelationSpearman"){ylim=c(0.55, 0.65)
  }else if (var == paste0("Tau", sp_name)){ylim=c(0.62, 0.75)
  }else if (var == "MeanRPKM"){ylim=c(6, 12)}
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

  
  for(enh in enhancer.datasets[[sp]]){
    for(class in levels(regland[["original"]][[enh]]$class_nb_contact)){
      this.genes=rownames(regland[["original"]][[enh]][which(regland[["original"]][[enh]]$class_nb_contact == class),])
      
      xpos=seq(1, length(levels(regland[["original"]][[enh]]$class_nb_contact)), 1)
      names(xpos) = levels(regland[["original"]][[enh]]$class_nb_contact)
      x=xpos[class]+smallx[enh]
      
      b=boxplot(expdiv[this.genes, var], plot=FALSE)
      med=median(expdiv[this.genes, var])
      ci=as.numeric(b$conf)
      
      points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, ci[1], x, ci[2],  col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
  mtext(names_MeasuresCM[plot.nb], side=2, line=2.5, cex=0.9)
  mtext(letters[plot.nb], side=3, line=1, at=0.1, font=2, cex=1.2)
  
  axis(side=1, cex.axis=1.2); mtext("Quantile of Complexity", side=1, line=2.5, cex=1)
  if (plot.nb == "1"){legend("topright", col=col.enhancers, legend = label.enhancers, box.col="white", bg="white", pch=20, cex=1, inset=c(-0.02,-0.05))}
}

######## Output ######

MeasuresCM <- c(paste0("Tau", sp_name), "MeanRPKM", "EuclideanSimilarity", "CorrelationSpearman",  "CorrectedEuclideanSimilarity", "CorrectedSpearman")
names_MeasuresCM <- c("Specificity (Tau)", "Expression level (RPKM)", "Euclidean Similarity", "Spearman's rho",  "Residual Euclidean Similarity", "Residual Spearman's rho" )

pdf(paste(pathFigures, "/SupplementaryFigure31.pdf", sep=""), width=7, height=7)
par(mfrow=c(3,2))
par(mai = c(0.5, 0.6, 0.3, 0.4)) # bottom, left, top, right

for (measure in 1:length(MeasuresCM)){
  CMPlot(MeasuresCM[measure], measure)
}

dev.off()

#############################################################################################################
######################## Gene expression level in Common cell types #########################################
CellTypesPlot <- function(var, plot.nb){
  if (var == "Conservation"){ylim=c(0.2, 0.6)}
  else if (var == "MeanRPKM"){ylim=c(6, 35)}
  else if (var == "ResidualConservation"){ylim=c(-0.15, 0.1)}
  
  for (cell in cells){
    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
    
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
    
    abline(v=xpos[1:4]+0.5, lty=3, col="gray40") # verticale lines between quantile
    
    if (plot.nb == "1"){mtext(cell, side=3, cex=1.2)} # Main in first row
    if (cell == "ESC" & plot.nb == "1"){legend("top", col=col.enhancers, legend = label.enhancers,pch=20, cex=1, bg="white", box.col="white")}
    
    if (cell == "Bcell"){axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1) # Yaxis labels in first column
        mtext(names_MeasuresCellTypes[plot.nb], side=2, line=2.5, cex=0.9)
        mtext(letters[plot.nb], side=3, line=1, at=0.1, font=2, cex=1.2)} 
      
        axis(side=1, cex.axis=1.2) # X axis labels in last row
        mtext("Complexity", side=1, line=2, cex=0.9)
    
  }
  
}

### output ###
MeasuresCellTypes <- c("MeanRPKM", "Conservation", "ResidualConservation") 
names_MeasuresCellTypes <- c("Expression level (RPKM)", "Expression Conservation", "Residual Expression Conservation")

pdf(paste(pathFigures, "/SupplementaryFigure30.pdf", sep=""), width=7, height=7)

par(mfrow=c(3,3))
par(mai = c(0.5, 0.5, 0.3, 0)) # bottom, left, top, right

for (measure in 1:length(MeasuresCellTypes)){
  CellTypesPlot(MeasuresCellTypes[measure], measure)
}

dev.off()