######################################################################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  library(Hmisc)
  source("../main_figures/parameters.R")
}

##############################################################################
if(load){
  sp="human"
  
  load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.RData", sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".regland.conservation.RData", sep=""))

  if (sp == "human"){
    sp_name="Human"
  } else{
    sp_name="Mouse"
  }

  expdiv$EuclideanSimilarity = 1-expdiv$EuclideanDistance
  distances =  "all" 

  load=FALSE
}


################################################################################################################################
############################## Cardoso-Moreira  - Euclidean Similarity ##########################################################

plot_profiles <- function(class_conserv, distances, xlab, xnames, plot.labels){
  smallx=c(-0.15, -0.075, 0.075, 0.15)
  names(smallx)=enhancer.datasets[[sp]]
  
  if (Measure == "corrected"){
    DivergenceMeasure = "CorrectedEuclideanSimilarity"
    ylab="1-Euclidean distance\n(corrected)"
    ylim=c(0.002, 0.02)
  }else{
    DivergenceMeasure = "EuclideanSimilarity"
    ylab="1-Euclidean distance"
    ylim=c(0.89, 0.94)
  }
  
  if (class_conserv == "class_cons_synt"){
    xmax=3
  }else{
    xmax=5
  }
  
  xlim=c(0.5, xmax+0.5)
  
  for(i in 1:length(distances)){

    dist=distances[i]
    
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
    
    axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
    mtext(ylab, side=2, line=2.5, cex=0.9)

     if (Measure == "corrected"){
       labelpos=xlim[1]-diff(xlim)/6
     } else{
       labelpos=xlim[1]-diff(xlim)/7
     }
    
    mtext(plot.labels[i], side=3, at=labelpos, line=1, font=2, cex=1.25)
    
  }
  
}

################################################################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure29.pdf", sep=""), width=6.85, height=10)
m=matrix(rep(NA, 4*2), nrow=4)

m[,1]=c(1:4)
m[,2]=c(5:8)

layout(m)

par(mar=c(4.1, 4.1, 2.1, 1)) # bottom, left, top, right

################################################################################################################################

## Gene expression profiles uncorrected
Measure = "uncorrected"

plot_profiles("class_nb_contact", distances,  "number of contacts class", 1:5, "a")

legend("bottomright", legend=enhancer.datasets[[sp]], pch=20,
       col=col.enhancers, cex=1, bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))

plot_profiles("class_align_score", distances,  "enhancer sequence conservation", 1:5, "b")
plot_profiles("class_cons_synt", distances,  "% conserved synteny", c("<75%", "75-99%", ">99%"), "c")
plot_profiles("class_cons_cont", distances,  "% conserved contacts", c("<1%", "25%", "50%", "75%", ">75%"), "d")

## Gene expression profiles corrected
Measure = "corrected"
plot_profiles("class_nb_contact", distances,  "number of contacts class", 1:5, "e")

plot_profiles("class_align_score", distances,  "enhancer sequence conservation", 1:5, "f")
plot_profiles("class_cons_synt", distances,  "% conserved synteny", c("<75%", "75-99%", ">99%"), "g")
plot_profiles("class_cons_cont", distances,  "% conserved contacts", c("<1%", "25%", "50%", "75%", ">75%"), "h")

dev.off()

################################################################################################################################
