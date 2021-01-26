###########################################################################################################

source("../main_figures/parameters.R") ## pathFinalData are defined based on the user name

sp="human"
sp_name="Human"

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/data.common.cells.expdiv.Rdata", sep=""))
load(paste(pathFigures, "RData/data.", sp, ".common.cells.regland.conservation.RData", sep=""))

if (sp == "human"){
  sp_name="Human"
}else{
  sp_name="Mouse"
}

cells <- c("Bcell", "ESC", "adipo")
col.cells = c("navy", "forestgreen", "darkorange")
names(col.cells) = cells
enh = "ENCODE"

#############################################################################################################
######################## Gene expression level in common cell types #########################################

xpos=seq(1, 5, 1)
smallx=c(-0.15, 0, 0.15)
names(smallx)=cells
xlim=c(0.5, 5.5)

CellTypesPlot <- function(var, plot.nb){
  if (var == "ExpressionConservation"){
    ylim=c(0.2, 0.6)
  }
  else if (var == paste0(sp, "_MeanRPKM")){
    ylim=c(0, 22)
  }
  else if (var == "ResidualExpressionConservation"){
    ylim=c(-0.11, 0.03)
  }
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

  for (cell in cells){
    regland = genes.conservation.cells[[enh]][[cell]][["obs"]][["all"]]
      genes = intersect(rownames(regland), rownames(expdiv_cells))
      regland = regland[genes,]
      
      for(class in levels(regland$class_nb_contact)){
        this.genes=rownames(regland[which(regland$class_nb_contact == class),])
        
        xpos=seq(1, length(levels(regland$class_nb_contact)), 1)
        names(xpos) = levels(regland$class_nb_contact)
        x=xpos[class]+smallx[cell]
        
        b=boxplot(expdiv_cells[this.genes, paste(cell, var, sep="_")], plot=FALSE)
        med=median(expdiv_cells[this.genes, paste(cell, var, sep="_")])
        ci=as.numeric(b$conf)
        
        points(x, med, pch=20, col=col.cells[cell], cex=1.1)
        segments(x, ci[1], x, ci[2],  col=col.cells[cell])
        
      }
    }
  
  abline(v=xpos[1:4]+0.5, lty=3, col="gray40") # verticale lines between quantile
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1) # Yaxis labels in first column
  mtext(names_MeasuresCellTypes[plot.nb], side=2, line=2.5, cex=0.75)
  
  ## plot labels
  mtext(letters[plot.nb], side=3, line=1, at=-0.8, font=2, cex=1.1)
  
  axis(side=1, cex.axis=1.1, mgp=c(3, 0.75, 0), at=1:5, labels=rep("",5)) # X axis labels in last row
  mtext(c("low", "medium", "high"), at=xpos[c(1,3,5)], side=1, line=0.75, cex=0.75)
  mtext("number of contacts", side=1, line=2, cex=0.75)
  
  if (plot.nb == 1){
    legend("topleft", legend=c("pre-adipocytes", "ESC", "B cells"), pch=20, col=rev(col.cells), cex=1.1, bty="o", 
           box.col="white", bg="white",  inset=c(0.01, -0.1), xpd=NA)
  }
}


################################################################################################################################
################### Gene expression evolution to Regulatory Landscape Evolution ################################################
n = 1

plot_cell <- function(class_conserv, Measure, distances, xlab, xnames){
  smallx=c(-0.15, 0, 0.15)
  names(smallx)=cells
  
  if (Measure == "corrected"){
    DivergenceMeasure = "ResidualExpressionConservation"
    ylab="expression conservation\n(corrected)"
    YLIM=c(-0.15, 0.11)
  } else{
    DivergenceMeasure = "ExpressionConservation"
    ylab="expression conservation"
    YLIM=c(0.1,0.75)
  } 
  
  if (class_conserv == "class_cons_synt"){
    xmax=3
  }else{
    xmax=5
  }
  
  xlim=c(0.5, xmax+0.5)
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=YLIM, xaxs="i", yaxs="i")
  
  for (cell in cells){
    regland = genes.conservation.cells[[enh]][[cell]][["obs"]][["all"]]
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

  if(class_conserv=="class_align_score"){
    mtext(c("low", "medium", "high"), at=xpos[c(1,3,5)], side=1, line=0.75, cex=0.75)
  } else{
    mtext(xnames, at=xpos, side=1, line=0.75, cex=0.75)
  }
  mtext(xlab, side=1, line=2.25, cex=0.75)
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1)
  mtext(ylab, side=2, line=2.5, cex=0.75)
  
}

################################################################################################################################
################################################### Output #####################################################################
# Relation with complexity
MeasuresCellTypes <- c(paste0(sp, "_MeanRPKM"), "ExpressionConservation", "ResidualExpressionConservation") 
names_MeasuresCellTypes <- c("expression level (TPM)", "expression conservation", "expression conservation\n(corrected)")

pdf(paste(pathFigures, "/ExtendedFigure3.pdf", sep=""), width=6.85, height=5.5)

par(mfrow=c(2,3))
par(mar = c(3.5, 4.65, 2.5, 1.1)) # bottom, left, top, right

for (measure in 1:length(MeasuresCellTypes)){
  CellTypesPlot(MeasuresCellTypes[measure], measure)
}

# Relation with Regulatory Landscape Evolution
distances =  "all"  # c("25kb - 100kb", "100kb - 500kb", "500kb - 2Mb", "all")

for (measure in c("corrected")){
  plot_cell("class_align_score", measure, distances, "enhancer sequence conservation", 1:5)
  mtext("d", side=3, at=-0.65, font=2, cex=1.1, line=1)

  plot_cell("class_cons_synt",  measure, distances, "% conserved synteny", c("<75%", "75-99%", ">99%"))
  mtext("e", side=3, at=-0.2, font=2, cex=1.1, line=1)
  
  plot_cell("class_cons_cont",  measure, distances, "% conserved contacts", c("<1%", "25%", "50%", "75", ">75%"))
  mtext("f", side=3, at=-0.65, font=2, cex=1.1, line=1)
} 

dev.off()

################################################################################################################################
