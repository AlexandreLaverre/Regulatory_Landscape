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
Measure = "corrected"
if (Measure == "corrected"){pdf_name="SupplementaryFigure33.pdf"}else{pdf_name="Figure6.pdf"}

#pdf(file=paste(pathFigures, pdf_name, sep=""), width = 8.5)

par(mfrow=c(1,4))
par(mai = c(0.5, 0.5, 0.5, 0)) # bottom, left, top, right

################################################################################################################################
######################################### PART1 : Common cells types  ##########################################################
if (Measure == "corrected"){
  DivergenceMeasure = "ResidualExpressionConservation"
  ylab="Residual Expression Level Conservation"
  YLIM=c(-0.2, 0.2)
}else{
  DivergenceMeasure = "ExpressionConservation"
  ylab="Expression Level Conservation"
  YLIM=c(0.35,0.85)
  } 

colnames(expdiv_cells)[3:5] <- c("adipo_ExpressionConservation", "Bcell_ExpressionConservation", "ESC_ExpressionConservation")


smallx=c(-0.15, 0, 0.15)
names(smallx)=cells

col.cells = c("navy", "forestgreen", "darkorange")
names(col.cells) = cells
enh = "ENCODE"
########################  A - Expression Conservation vs Number of conserved enhancers ######################## 

plot_cell <- function(class_conserv, xlab, xnames){
  par(mfrow=c(1,4))
  #par(mai = c(0.5, 0.5, 0.5, 0)) # bottom, left, top, right
  if (class_conserv == "class_cons_synt"){xmax=2}else{xmax=5}
  
  xlim=c(0.5, xmax+0.5)
  
  for (dist in names(genes.conservation.cells[[enh]][["ESC"]][["obs"]])){
    
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
    
    if (dist == "25kb - 100kb"){
      axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
      mtext(ylab, side=2, line=2.5, cex=0.9)
      mtext("d", side=3, at=0.45, font=2, cex=1.2, line=0.5)
      legend("bottomleft", legend=cells, pch=20,
             col=col.cells, cex=1.2, bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))
    }
    
    mtext(dist, side=3, line=-1, cex=1)
  }
  
  ## plot label
  mtext("a", side=3, line=1, at=1, font=2, cex=1.05)
}

plot_cell("class_align_score", "Alignement score", 1:5)
plot_cell("class_cons_synt", "Synteny Conservation", c("without", "with"))
plot_cell("class_cons_cont", "Contact conservation", c("<1%", "25%", "50%", "75", ">75%"))

# ################################################################################################################################
# ############################## PART2 : All cells & Cardoso-Moreira  ##########################################################
# pdf(paste(pathFigures, "/Fig6_by_distances_all.pdf", sep=""), width = 7, height = 3)
# 
# par(mfrow=c(1,3))
# par(mai = c(0.5, 0.5, 0.2, 0.2)) # bottom, left, top, right
# Measure = "uncorrected"
# data = "obs"
# 
# if (Measure == "corrected"){DivergenceMeasure = "CorrectedSpearman"; xlab="Residual Spearman's rho"
# }else{DivergenceMeasure = "CorrelationSpearman"; xlab="Spearman's rho"} 
# 
# #if (DivergenceMeasure == "CorrectedEuclideanSimilarity"){ylim=c(0, 0.02)}else{ylim=c(0.005, 0.045)} #-0.05, 0.15
# if (DivergenceMeasure == "CorrectedSpearman"){ylim=c(-0.02, 0.15)}else{ylim=c(0.5, 0.68)} #
# 
# #### D - Gene expression profil similarity and enhancers conserved in sequences ####
# first_enh="ENCODE"
# first_dist = names(genes.conservation[[first_enh]][[data]])[1]
# 
# smallx=c(-0.15, -0.075, 0.075, 0.15)
# names(smallx)=enhancer.datasets[[sp]]
# 
# xlim=c(0.5, length(levels(genes.conservation[[first_enh]][[data]][[first_dist]]$class_align_score))+0.5)
# 
# for (dist in names(genes.conservation[[first_enh]][[data]])){
#   
#   plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
#   
#   for(enh in enhancer.datasets[[sp]]){
#     message(enh)
#     regland = genes.conservation[[enh]][[data]][[dist]]
#     genes = intersect(rownames(regland), rownames(expdiv))
#     regland = regland[genes,]
#     
#     for(class in levels(regland$class_align_score)){
#       this.genes=rownames(regland[which(regland$class_align_score == class),])
#       
#       xpos=seq(1, length(levels(regland$class_align_score)), 1)
#       names(xpos) = levels(regland$class_align_score)
#       x=xpos[class]+smallx[enh]
#       
#       b=boxplot(expdiv[this.genes, DivergenceMeasure], plot=FALSE)
#       med=median(expdiv[this.genes, DivergenceMeasure])
#       ci=as.numeric(b$conf)
#       
#       points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
#       segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
#     }
#   }
#   
#   abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
#   axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland$class_align_score))), cex.axis=0.8)
#   mtext(xpos, at=xpos, side=1, line=1, cex=0.8)
#   mtext("Alignment score", side=1, line=2.5, cex=0.9)
#   
#   if (dist == first_dist){
#     axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
#     mtext(xlab, side=2, line=2.5, cex=0.9)
#   }
#   
#   # if (dist == first_dist){
#   #   mtext("d", side=3, at=0.45, font=2, cex=1.05, line=0.5)
#   #   legend("bottomleft", legend=label.enhancers[enhancer.datasets[[sp]]], pch=20,
#   #          col=col.enhancers[enhancer.datasets[[sp]]], cex=1,
#   #          bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))
#   # }
#   
#   mtext(dist, side=3, line=-1, cex=1)
# }
# 
# 
# #### E - Gene expression profil similarity and enhancers conserved in synteny ####
# if (DivergenceMeasure == "CorrectedSpearman"){ylim=c(-0.02, 0.15)}else{ylim=c(0.5, 0.7)} #
# 
# 
# xlim=c(0.5, 2.5)
# 
# for (dist in names(genes.conservation[[enh]][[data]])){
#   
#   plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
#   
#   for(enh in enhancer.datasets[[sp]]){
#     message(enh)
#     regland = genes.conservation[[enh]][[data]][[dist]]
#     genes = intersect(rownames(regland), rownames(expdiv))
#     regland = regland[genes,]
#     
#     for(class in levels(regland$class_cons_synt)){
#       this.genes=rownames(regland[which(regland$class_cons_synt == class),])
#       
#       xpos=seq(1, length(levels(regland$class_cons_synt)), 1)
#       names(xpos) = levels(regland$class_cons_synt)
#       x=xpos[class]+smallx[enh]
#       
#       b=boxplot(expdiv[this.genes, DivergenceMeasure], plot=FALSE)
#       med=median(expdiv[this.genes, DivergenceMeasure])
#       ci=as.numeric(b$conf)
#       
#       points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
#       segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
#     }
#   }
#   
#   abline(v=xpos[1]+0.5, lty=3, col="gray40")
#   axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland$class_cons_synt))), cex.axis=0.8)
#   mtext(c("<99%", ">99%"), at=xpos, side=1, line=1, cex=0.8)
#   #mtext("Conservation of synteny", side=1, line=2.5, cex=0.9)
#   
#   if (dist == first_dist){
#     axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
#     mtext(xlab, side=2, line=2.5, cex=0.9)
#   }
#   
#   # if (dist == first_dist){
#   #   mtext("e", side=3, at=0.45, font=2, cex=1.05, line=0.5)
#   #   legend("bottomleft", legend=label.enhancers[enhancer.datasets[[sp]]], pch=20,
#   #          col=col.enhancers[enhancer.datasets[[sp]]], cex=1,
#   #          bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))
#   # }
#   # 
#   # mtext(dist, side=3, line=-5, cex=1.2)
# }
# 
# #### F - Gene expression profil similarity and enhancers conserved in contact ####
# if (DivergenceMeasure == "CorrectedSpearman"){ylim=c(-0.02, 0.15)}else{ylim=c(0.5, 0.7)} #
# 
# xlim=c(0.5, 5.5)
# 
# for (dist in names(genes.conservation[[enh]][[data]])){
#   plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
#   
#   for(enh in enhancer.datasets[[sp]]){
#     message(enh)
#     regland = genes.conservation[[enh]][[data]][[dist]]
#     genes = intersect(rownames(regland), rownames(expdiv))
#     regland = regland[genes,]
#     
#     for(class in levels(regland$class_cons_cont)){
#       this.genes=rownames(regland[which(regland$class_cons_cont == class),])
#       
#       xpos=seq(1, length(levels(regland$class_cons_cont)), 1)
#       names(xpos) = levels(regland$class_cons_cont)
#       x=xpos[class]+smallx[enh]
#       
#       b=boxplot(expdiv[this.genes, DivergenceMeasure], plot=FALSE)
#       med=median(expdiv[this.genes, DivergenceMeasure])
#       ci=as.numeric(b$conf)
#       
#       points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
#       segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
#     }
#   }
#   
#   abline(v=xpos[1:3]+0.5, lty=3, col="gray40")
#   axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland$class_cons_cont))), cex.axis=0.8)
#   mtext(c("<1%", "25%", "50%", "75", ">75%"), at=xpos, side=1, line=1, cex=0.7)
#   mtext("Conservation of contact", side=1, line=2.5, cex=0.9)
#   
#   if (dist == first_dist){
#     axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
#     mtext(xlab, side=2, line=2.5, cex=0.9)
#   }
#   
#   # if (dist == first_dist){
#   #   mtext("e", side=3, at=0.45, font=2, cex=1.05, line=0.5)
#   #   legend("bottomleft", legend=label.enhancers[enhancer.datasets[[sp]]], pch=20,
#   #          col=col.enhancers[enhancer.datasets[[sp]]], cex=1,
#   #          bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))
#   # }
#   
#   # mtext(dist, side=3, line=-5, cex=1.2)
# }
# 
# dev.off()