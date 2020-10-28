#########################################################################################################################
source("parameters.R") ## pathFiguress are defined based on the user name

library(ape)
library(vioplot)

ref_sp = "human"
target_sp = "mouse"

load(paste(pathFigures, "RData/Fig4_", ref_sp, ".Rdata", sep=""))

species <- c("macaque", "dog", "cow", "elephant", "rabbit", "rat", target_sp, "opossum", "chicken")

#########################################################################################################################
if(ref_sp == "human"){pdf_name="Figure4.pdf"}else{pdf_name="Sup_Figure12.pdf"}

pdf(paste(pathFigures, pdf_name, sep=""))
par(mai = c(0.7, 0.8, 0.4, 0.2)) # bottom, left, top, right
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE))

############################################   A - Global synteny conservation ############################################ 
# Gene - ENCODE synteny
enh="ENCODE"
bar <- barplot((1-conserv_synteny[[enh]]$data)*100, beside = TRUE, ylim = c(0,30), xpd=FALSE, las=2,
               border=c(dataset.colors, "white"), col=c(dataset.colors, "white"),
               lwd=1.5, cex.names=0.8, density=c(dataset.density,0), angle=c(dataset.angle,0),
               main="", ylab='Re-arranged \n gene-enhancer pairs (%)')

# Axis and legend
box(bty="l")
bar_position = seq(1.3,33,3.6)
middle_bar = seq(2,27,3)

axis(side=1, at=bar_position, labels=species, mgp=c(3, 0.65, 0), cex.axis=1.1, las=2)

legend("topleft", legend = c("Original", "Simulated"), 
       border=dataset.colors, density=c(dataset.density,0), angle=c(dataset.angle,0), bty='n', cex=0.8)

# Confidence intervals
arrows(x0=bar,y0=(1-conserv_synteny[[enh]]$conf_up)*100,y1=(1-conserv_synteny[[enh]]$conf_low)*100,angle=90,code=3,length=0.05)
segments(x0=bar_position-0.6,x1=bar_position+0.6, y0=(1-conserv_synteny[[enh]]$data[middle_bar]+0.02)*100)

# Significance of the test
x=1
for (test in conserv_synteny[[enh]]$p_test){
  if (test < 0.0001){text("***",x = bar_position[x], y=(1-conserv_synteny[[enh]]$data[middle_bar[x]]+0.03)*100, cex=1)}
  else if (test < 0.001){text("**",x = bar_position[x], y=(1-conserv_synteny[[enh]]$data[middle_bar[x]]+0.03)*100, cex=1)}
  else if (test < 0.01){text("*",x = bar_position[x], y=(1-conserv_synteny[[enh]]$data[middle_bar[x]]+0.03)*100, cex=1)}
  else {text("NS",x = bar_position[x], y=(1-conserv_synteny[[enh]]$data[middle_bar[x]]+0.03)*100, cex=0.8)}
  x = x + 1
}

mtext("A", side=3, line=1, at=-1.5, font=2, cex=1.2)

####################################  B & C - Synteny conservation and distance from promoters - all enhancers datasets ############################### 
species = c(target_sp, "chicken")
par(mai = c(0.7, 0.8, 0.4, 0.1)) # bottom, left, top, right

for (sp in species){
  for (enh in enhancer.datasets[[ref_sp]]){
    conserv = get(paste("conserv_synteny_dist_", enh, sep=""))
    
    if (enh == "ENCODE"){ # First sp
      if (sp == target_sp){
        YLAB="Excess of gene-enh \n conserved in synteny (%)"
        YLIM=c(-2,15)
      }else{YLAB=""
      YLIM=c(-10,80)
      par(mai = c(0.7, 0.4, 0.4, 0.4))} # bottom, left, top, right
      
      plot((conserv[[sp]]$obs-conserv[[sp]]$simul)*100/conserv[[sp]]$simul, type="l",
           col=col.enhancers[[enh]], xaxt = "n", las=2, lwd=0.8,
           xlab="", ylab=YLAB, ylim=YLIM, main="")
      
      # Add lines of other enhancers datasets
    }else{lines((conserv[[sp]]$obs-conserv[[sp]]$simul)*100/conserv[[sp]]$simul, type="l", col=col.enhancers[[enh]] , lwd=0.8)} 
    
    # Add confidences intervals
    for (row in 1:nrow(conserv[[sp]])){
      segments(x0=row,y0=(conserv[[sp]][row,"conf_low_obs"]-conserv[[sp]][row,"conf_low_simul"])*100/conserv[[sp]][row,"conf_low_simul"],
               x1=row,y1=(conserv[[sp]][row,"conf_up_obs"]-conserv[[sp]][row,"conf_up_simul"])*100/conserv[[sp]][row,"conf_up_simul"],
               col=col.enhancers[[enh]], lwd=0.3)
    }
    
  }
  
  # Axis and legends
  abline(h=0, lty=2, col='black')
  axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
  mtext("Distance from TSS (Mb)", side=1, line=2.25, cex=0.85)
  
  if (sp==target_sp){mtext("B", side=3, line=1, at=-1.5, font=2, cex=1.2)
    legend("topleft", col=col, legend = enhancers, bty='n', lty=1, cex=0.8)}else{mtext("C", side=3, line=1, at=-1.5, font=2, cex=1.2)}
  
  
}


dev.off()
