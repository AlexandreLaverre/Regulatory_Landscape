#########################################################################################################################
source("parameters.R") ## pathFiguress are defined based on the user name

library(ape)
library(vioplot)

ref_sp = "human"
target_sp = "mouse"

load(paste(pathFigures, "RData/Fig4_", ref_sp, ".Rdata", sep=""))

species <- c("macaque", "dog", "cow", "elephant", "rabbit", "rat", target_sp, "opossum", "chicken")

enhancers <- c("FANTOM5", "ENCODE")
if(ref_sp == "human"){enhancers <- c(enhancers, "RoadMap", "GRO_seq")}

#########################################################################################################################
if(ref_sp == "human"){pdf_name="Figure4.pdf"}else{pdf_name="Sup_Figure12.pdf"}

#pdf(paste(pathFigures, "Figure4_", ref_sp, ".pdf", sep=""), width=8.5, height=8)
par(mai = c(0.7, 0.8, 0.4, 0.2)) # bottom, left, top, right
layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))

############################################   A - Global synteny conservation ############################################ 
# Gene - ENCODE synteny
enh="ENCODE"
bar <- barplot(conserv_synteny[[enh]]$data, beside = TRUE, ylim = c(0.6,1.03), xpd=FALSE, las=2,
               border = c("darkgreen", "firebrick3", "white"),col="white",
               main="", ylab='Proportion of gene-enhancer maintained in synteny')

# Axis and legend
box(bty="l")
bar_position = seq(1.3,33,3.6)
middle_bar = seq(2,27,3)

axis(side=1, at=bar_position, labels=species, mgp=c(3, 0.65, 0), cex.axis=1.1, las=2)

legend("topright", border=c("darkgreen", "firebrick3"), fill="white", legend = c("Original", "Simulated"), bty='n', cex=0.8)

# Confidence intervals
arrows(x0=bar,y0=conserv_synteny[[enh]]$conf_up,y1=conserv_synteny[[enh]]$conf_low,angle=90,code=3,length=0.05)
segments(x0=bar_position-0.6,x1=bar_position+0.6, y0=conserv_synteny[[enh]]$data[middle_bar]+0.02)

# Significance of the test
x=1
for (test in conserv_synteny[[enh]]$p_test){
  if (test < 0.0001){text("***",x = bar_position[x], y=conserv_synteny[[enh]]$data[middle_bar[x]]+0.03, cex=1)}
  else if (test < 0.001){text("**",x = bar_position[x], y=conserv_synteny[[enh]]$data[middle_bar[x]]+0.03, cex=1)}
  else if (test < 0.01){text("*",x = bar_position[x], y=conserv_synteny[[enh]]$data[middle_bar[x]]+0.03, cex=1)}
  else {text("NS",x = bar_position[x], y=conserv_synteny[[enh]]$data[middle_bar[x]]+0.03, cex=0.8)}
  x = x + 1
}

mtext("A", side=3, line=1, at=-1.5, font=2, cex=1.2)

###################################  B - Synteny conservation and distance from promoters - all species ################################# 
if (ref_sp == "human"){YMAX=0.25}else{YMAX=0.35}

conserv = get(paste("conserv_synteny_dist", enh, sep="_"))
col <- c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#ABD9E9","#74ADD1","#4575B4","#313695")
color_n = 1 # To change color between each enhancers dataset

par(lwd = 1.2) 
for (sp in species){
  # First sp
  if (sp == "macaque"){ 
    plot((conserv[[sp]]$obs-conserv[[sp]]$simul)/conserv[[sp]]$simul, type="l", col=col[color_n], xaxt = "n", las=2, ylim=c(-0.2,YMAX), 
         xlab="", ylab="Original-Simulated \n gene-enh maintained in synteny", main="", lwd=0.8)
  
  # Add lines of other species
  }else{lines((conserv[[sp]]$obs-conserv[[sp]]$simul)/conserv[[sp]]$simul, type="l", col=col[color_n], lwd=0.8)} 
  
  # Add confidences intervals
  for (row in 1:nrow(conserv[[sp]])){
    segments(x0=row,y0=(conserv[[sp]][row,"conf_low_obs"]-conserv[[sp]][row,"conf_low_simul"])/conserv[[sp]][row,"conf_low_simul"],
             x1=row,y1=(conserv[[sp]][row,"conf_up_obs"]-conserv[[sp]][row,"conf_up_simul"])/conserv[[sp]][row,"conf_up_simul"], col=col[color_n], lwd=0.3)}
  
  color_n = color_n + 1 
}

# Axis and legends
par(lwd = 1) 
abline(h=0, lty=2, col='black')
class_leg <- c("0", "0.5", "1", "1.5", "2")
axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
mtext("Distance from TSS (Mb)", side=1, line=2.25, cex=0.85)
legend("bottomleft", col=col, ncol=3, legend = species, bty='n', lty=1, cex=0.8)
mtext("B", side=3, line=1, at=-1.5, font=2, cex=1.2)

####################################  C & D - Synteny conservation and distance from promoters - all enhancers datasets ############################### 
species = c(target_sp, "dog")

for (sp in species){
  col <- c("red", "navy", "forestgreen", "orange")
  color_n = 1 # To change color between each enhancers dataset
  
  for (enhancer in enhancers){
    conserv = get(paste("conserv_synteny_dist_", enhancer, sep=""))
    
    if (enhancer == "FANTOM5"){ # First sp
      plot((conserv[[sp]]$obs-conserv[[sp]]$simul)/conserv[[sp]]$simul, type="l", col=col[color_n], xaxt = "n", las=2, lwd=0.8,
           xlab="", ylab="Original-Simulated \n gene-enh maintained in synteny", main="")
      
      # Add lines of other enhancers datasets
    }else{lines((conserv[[sp]]$obs-conserv[[sp]]$simul)/conserv[[sp]]$simul, type="l", col=col[color_n] , lwd=0.8)} 
    
    # Add confidences intervals
    for (row in 1:nrow(conserv[[sp]])){
      segments(x0=row,y0=(conserv[[sp]][row,"conf_low_obs"]-conserv[[sp]][row,"conf_low_simul"])/conserv[[sp]][row,"conf_low_simul"],
               x1=row,y1=(conserv[[sp]][row,"conf_up_obs"]-conserv[[sp]][row,"conf_up_simul"])/conserv[[sp]][row,"conf_up_simul"],
               col=col[color_n], lwd=0.3)
    }
    
    color_n = color_n + 1 
  }
  
  # Axis and legends
  abline(h=0, lty=2, col='black')
  axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
  mtext("Distance from TSS (Mb)", side=1, line=2.25, cex=0.85)
  legend("topleft", col=col, legend = enhancers, bty='n', lty=1, cex=0.8)
  
  if (sp==target_sp){mtext("C", side=3, line=1, at=-1.5, font=2, cex=1.2)}else{mtext("D", side=3, line=1, at=-1.5, font=2, cex=1.2)}
  
  
}


#dev.off()
