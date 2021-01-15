###########################################################################################################

source("../main_figures/parameters.R") ## pathFinalData are defined based on the user name

sp="human"
sp_name="Human"

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/data.", sp, ".common.cells.expdiv.Rdata", sep=""))
load(paste(pathFigures, "RData/data.", sp, ".common.cells.regland.conservation.RData", sep=""))

expdiv$EuclideanSimilarity <- 1-expdiv$EuclideanDistance

if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}

cells <- c("Bcell", "ESC", "adipo")
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

#############################################################################################################
######################## Gene expression level in Common cell types #########################################
CellTypesPlot <- function(var, plot.nb){
  if (var == "ExpressionConservation"){ylim=c(0.2, 0.6)}
  else if (var == paste0(sp, "_MeanRPKM")){ylim=c(0, 25)}
  else if (var == "ResidualExpressionConservation"){ylim=c(-0.12, 0.05)}
  
  for (cell in cells){
    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
    
    for(enh in enhancer.datasets[[sp]]){
      regland = genes.conservation.cells[[enh]][[cell]][["obs"]][["all"]]
      genes = intersect(rownames(regland), rownames(expdiv_cells))
      regland = regland[genes,]
      
      for(class in levels(regland$class_nb_contact)){
        this.genes=rownames(regland[which(regland$class_nb_contact == class),])
        
        xpos=seq(1, length(levels(regland$class_nb_contact)), 1)
        names(xpos) = levels(regland$class_nb_contact)
        x=xpos[class]+smallx[enh]
        
        b=boxplot(expdiv_cells[this.genes, paste(cell, var, sep="_")], plot=FALSE)
        med=median(expdiv_cells[this.genes, paste(cell, var, sep="_")])
        ci=as.numeric(b$conf)
        
        points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
        segments(x, ci[1], x, ci[2],  col=col.enhancers[enh])
        
      }
    }
    
    abline(v=xpos[1:4]+0.5, lty=3, col="gray40") # verticale lines between quantile
    
    if (plot.nb == "1"){mtext(cell, side=3, cex=1.2)} # Main in first row
    if (cell == "Bcell" & plot.nb == "2"){legend("topleft", col=col.enhancers, legend = label.enhancers,pch=20, cex=1, bg="white", box.col="white")}
    
    if (cell == "Bcell"){axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1) # Yaxis labels in first column
      mtext(names_MeasuresCellTypes[plot.nb], side=2, line=2.5, cex=0.9)
      mtext(letters[plot.nb], side=3, line=1, at=0.1, font=2, cex=1.2)} 
    
    axis(side=1, cex.axis=1.2) # X axis labels in last row
    mtext("Quantile of Complexity", side=1, line=2, cex=0.9)
    
  }
  
}

### output ###
MeasuresCellTypes <- c(paste0(sp, "_MeanRPKM"), "ExpressionConservation", "ResidualExpressionConservation") 
names_MeasuresCellTypes <- c("Expression level (RPKM)", "Expression Conservation", "Residual Expression Conservation")

pdf(paste(pathFigures, "/SupplementaryFigure29.pdf", sep=""), width=7, height=7)

par(mfrow=c(3,3))
par(mai = c(0.5, 0.5, 0.3, 0)) # bottom, left, top, right

for (measure in 1:length(MeasuresCellTypes)){
  CellTypesPlot(MeasuresCellTypes[measure], measure)
}

dev.off()

