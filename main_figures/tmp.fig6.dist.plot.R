######################################################################################################################
library(Hmisc)
setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")

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
  
  load(paste(pathFigures, "RData/data.", sp, ".common.cells.expdiv.Rdata", sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".gene.regland.conservation.common.cells.RData", sep=""))
  
  if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}
  
  cells <- c("Bcell", "ESC", "adipo")
  dataset.colors=c("firebrick1", "forestgreen", "navy")
  names(dataset.colors) = cells
  load=FALSE
  
}

################################################################################################################################
Measure = "uncorrected"
if (Measure == "corrected"){pdf_name="SupplementaryFigure33.pdf"}else{pdf_name="Figure6.pdf"}

#pdf(file=paste(pathFigures, pdf_name, sep=""), width = 8.5)

par(mfrow=c(1,3))
par(mai = c(0.5, 0.5, 0.5, 0)) # bottom, left, top, right

################################################################################################################################
######################################### PART1 : Common cells types  ##########################################################
if (Measure == "corrected"){DivergenceMeasure = "ResidualExpressionConservation"; ylab="Residual Expression Level Conservation";
YLIM=c(-0.15, 0.11)
}else{DivergenceMeasure = "ExpressionConservation"; ylab="Expression Level Conservation"; YLIM=c(0.1,0.75)} 


col.cells = c("navy", "forestgreen", "darkorange")
names(col.cells) = cells
enh = "ENCODE"

plot_cell <- function(class_conserv, xlab, xnames){
  smallx=c(-0.15, 0, 0.15)
  names(smallx)=cells
  
  #par(mai = c(0.5, 0.5, 0.5, 0)) # bottom, left, top, right
  if (class_conserv == "class_cons_synt"){xmax=2}else{xmax=5}
  xlim=c(0.5, xmax+0.5)
  
  for (dist in rev(names(genes.conservation.cells[[enh]][["ESC"]][["obs"]]))){
    
    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM, xaxs="i", yaxs="i")
    
    for (cell in cells){
      regland = genes.conservation.cells[[enh]][[cell]][["obs"]][[dist]]
      genes = intersect(rownames(regland), rownames(expdiv_cells))
      regland = regland[genes,]
      
      for (class in levels(regland[[class_conserv]])){
        this.genes=rownames(regland[which(regland[[class_conserv]] == class),])
        
        b=boxplot(expdiv_cells[this.genes, paste(cell, "_", DivergenceMeasure, sep="")], plot=F)
        med=median(expdiv_cells[this.genes, paste(cell, "_", DivergenceMeasure, sep="")])
        ci=as.numeric(b$conf)
        
        xpos=seq(1,  length(levels(regland[[class_conserv]])), 1)
        names(xpos) = levels(regland[[class_conserv]])
        
        x=xpos[class]+smallx[cell]
        
        points(x, med, pch=20, col=col.cells[cell], cex=1.1)
        segments(x, ci[1], x, ci[2], col=col.cells[cell])
      }
    }
    
    abline(v=xpos[1:xmax-1]+0.5, lty=3, col="gray40")
    axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[class_conserv]]))), cex.axis=0.8)
    
    mtext(xnames, at=xpos, side=1, line=1, cex=0.8)
    mtext(xlab, side=1, line=2.5, cex=0.9)
    
    if (dist == "all"){
      axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
      mtext(ylab, side=2, line=2.5, cex=0.9)
      mtext("d", side=3, at=0.45, font=2, cex=1.2, line=0.5)
      legend("topleft", legend=rev(cells), pch=20,
             col=rev(col.cells), cex=1.2, bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))
      
    }

    mtext(dist, side=3, line=-1, cex=1)
  }

}

par(mfrow=c(1,4))
plot_cell("class_align_score", "Alignement score", 1:5)
# mtext("d", side=3, at=0.45, font=2, cex=1.2, line=0.5)
# legend("topleft", legend=rev(cells), pch=20,
#        col=rev(col.cells), cex=1.2, bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))

plot_cell("class_cons_synt", "Synteny Conservation", c("without", "with"))
plot_cell("class_cons_cont", "Contact conservation", c("<1%", "25%", "50%", "75", ">75%"))



# ################################################################################################################################
# ############################## PART2 :  Cardoso-Moreira  ##########################################################
#pdf(paste(pathFigures, "/Fig6_by_distances_all.pdf", sep=""), width = 7, height = 3)

#Measure = "uncorrected"

if (Measure == "corrected"){DivergenceMeasure = "CorrectedSpearman"; ylab="Residual Spearman's rho"; ylim=c(-0.02, 0.13)
}else{DivergenceMeasure = "CorrelationSpearman"; ylab="Spearman's rho"; ylim=c(0.5, 0.68)}

# expdiv$EuclideanSimilarity = 1-expdiv$EuclideanDistance
# DivergenceMeasure = "CorrectedEuclideanSimilarity"
# 
# if (DivergenceMeasure == "EuclideanSimilarity"){ylab="Euclidean Similarity"; ylim=c(0.89, 0.94)
# }else{ylab="Residual Euclidean Similarity"; ylim=c(0, 0.02)}

plot_profiles <- function(class_conserv, xlab, xnames){
  #par(mfrow=c(1,4))
  #par(mai = c(0.5, 0.5, 0.5, 0)) # bottom, left, top, right
  
  smallx=c(-0.15, -0.075, 0.075, 0.15)
  names(smallx)=enhancer.datasets[[sp]]
  
  if (class_conserv == "class_cons_synt"){xmax=2}else{xmax=5}
  xlim=c(0.5, xmax+0.5)
  
  for (dist in rev(names(genes.conservation[[enh]][["obs"]]))[1]){
    
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
      mtext("d", side=3, at=0.45, font=2, cex=1.2, line=0.5)
      legend("bottomleft", legend=enhancer.datasets[[sp]], pch=20,
             col=col.enhancers, cex=1, bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))
    }
    
    mtext(dist, side=3, line=-1, cex=1)
  }
}

#par(mfrow=c(1,3))
plot_profiles("class_align_score", "Alignement score", 1:5)
plot_profiles("class_cons_synt", "Synteny Conservation", c("without", "with"))
plot_profiles("class_cons_cont", "Contact conservation", c("<1%", "25%", "50%", "75", ">75%"))


#dev.off()