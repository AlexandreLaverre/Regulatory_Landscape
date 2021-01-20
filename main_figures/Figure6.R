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
  load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.Rdata", sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".regland.conservation.RData", sep=""))
  
  if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}
  
  load=FALSE
  
}

#############################################################################################################
######################## Part 1 : Complexity and gene expression structure ##################################
xpos=seq(1, 5, 1)
smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]
xlim=c(0.5, 5.5)
cex.mtext = 0.8

CMPlot <- function(var, plot.nb){
  if (var == "CorrelationSpearman"){ylim=c(0.55, 0.65)
  }else if (var == paste0("Tau", sp_name)){ylim=c(0.60, 0.8)
  }else if (var == "MeanRPKM"){ylim=c(7, 11)
  }else if (var == "CorrectedSpearman"){ylim=c(0.02, 0.1)}
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    regland = genes.conservation[[enh]][["obs"]][["all"]]
    genes = intersect(rownames(regland), rownames(expdiv))
    regland = regland[genes,]
    
    for(class in levels(regland$class_nb_contact)){
      this.genes=rownames(regland[which(regland$class_nb_contact == class),])
      
      xpos=seq(1, length(levels(regland$class_nb_contact)), 1)
      names(xpos) = levels(regland$class_nb_contact)
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
  mtext(names_MeasuresCM[plot.nb], side=2, line=2.5, cex=cex.mtext)
  mtext(letters[plot.nb], side=3, line=1, at=0.1, font=2, cex=1.2)
  

  axis(side=1, cex.axis=1.2); mtext("Number of contacts (quantile)", side=1, line=2.5, cex=cex.mtext)
  if (plot.nb == "1"){legend("topleft", col=col.enhancers, legend = label.enhancers, box.col="white", bg="white",
                             pch=20, cex=0.9, inset=c(0.01, -0.18), xpd=T)}
}

################################################################################################################################
############################## PART2 :  Cardoso-Moreira plot function ##########################################################

plot_profiles <- function(class_conserv, distances, xlab, xnames){
  smallx=c(-0.15, -0.075, 0.075, 0.15)
  names(smallx)=enhancer.datasets[[sp]]
  
  if (Measure == "corrected"){DivergenceMeasure = "CorrectedSpearman"; ylab="Corrected Spearman's rho"; ylim=c(0, 0.13)
  }else{DivergenceMeasure = "CorrelationSpearman"; ylab="Spearman's rho"; ylim=c(0.5, 0.7)}
  
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
    axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[class_conserv]]))), cex.axis=cex.mtext)
    mtext(xnames, at=xpos, side=1, line=1, cex=cex.mtext)
    mtext(xlab, side=1, line=2.5, cex=cex.mtext)
    
    if (dist == "all"){
      axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
      mtext(ylab, side=2, line=2.5, cex=cex.mtext)
    }
    
  }
}

################################################################################################################################
################################################### Plot Figure 6 ##############################################################
Measure = "corrected"

pdf(file=paste(pathFigures, "Figure6.pdf", sep=""), width = 6.85,  heigh=5.5)

par(mfrow=c(2,3))
par(mai = c(0.5, 0.5, 0.5, 0.1))

######################## Part 1 : Complexity ######################## 
MeasuresCM <- c("MeanRPKM", paste0("Tau", sp_name),  "CorrectedSpearman")
names_MeasuresCM <- c("Mean Expression level (RPKM)", "Specificity", "Corrected Spearman's rho")

for (measure in 1:length(MeasuresCM)){
  CMPlot(MeasuresCM[measure], measure)
}

######################## Part 2 : Gene expression profiles evolution ################
plot_profiles("class_align_score", "all",  "Mean Alignement score (quantile)", 1:5)
mtext("d", side=3, at=0.45, font=2, cex=1.2, line=0.5)

plot_profiles("class_cons_synt", "all",  "Synteny Conservation", c("<75%", "75-99%", ">99%"))
plot_profiles("class_cons_cont", "all",  "Contact conservation", c("<1%", "25%", "50%", "75%", ">75%"))


dev.off()

