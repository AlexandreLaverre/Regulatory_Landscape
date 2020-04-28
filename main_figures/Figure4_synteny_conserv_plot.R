############################################################ PLOT FIGURE 5 ###############################################
path <- "/home/laverre/Data/Regulatory_landscape/result/Main_figures/"
load(paste(path, "Fig4_synteny_human.Rdata", sep=""))


pdf(paste(path, "Figure4_synteny_human.pdf", sep=""), width=8.5, height=8)
par(mai = c(0.8, 0.8, 0.3, 0.2)) # bottom, left, top, right
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))

enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")

############################################   A - Global synteny conservation ############################################ 
enh="ENCODE"

if (enh == "CAGE"){YMIN=0.6}else{YMIN=0.75}

bar <- barplot(conserv_synteny[[enh]]$data, beside = TRUE, ylim = c(YMIN,1.03), xpd=FALSE, las=2,
               border = c("darkgreen", "firebrick3", "white"),col="white",
               main="Fig.4.A", ylab='Proportion of gene-enhancer maintained in synteny')

box(bty="l")
bar_position = seq(1.3,33,3.6)
middle_bar = seq(2,27,3)

axis(1, at=bar_position, labels=F)
text(bar_position, par("usr")[3]-0.01, labels = species, pos = 1, xpd = TRUE, srt=20, cex=0.8)
legend("topright", border=c("darkgreen", "firebrick3"), fill="white", legend = c("Original", "Simulated"), bty='n', cex=0.8)

arrows(x0=bar,y0=conserv_synteny[[enh]]$conf_up,y1=conserv_synteny[[enh]]$conf_low,angle=90,code=3,length=0.05)
segments(x0=bar_position-0.6,x1=bar_position+0.6, y0=conserv_synteny[[enh]]$data[middle_bar]+0.02)

x=1
for (test in conserv_synteny[[enh]]$p_test){
  if (test < 0.0001){text("***",x = bar_position[x], y=conserv_synteny[[enh]]$data[middle_bar[x]]+0.03, cex=1)}
  else if (test < 0.001){text("**",x = bar_position[x], y=conserv_synteny[[enh]]$data[middle_bar[x]]+0.03, cex=1)}
  else if (test < 0.01){text("*",x = bar_position[x], y=conserv_synteny[[enh]]$data[middle_bar[x]]+0.03, cex=1)}
  else {text("NS",x = bar_position[x], y=conserv_synteny[[enh]]$data[middle_bar[x]]+0.03, cex=0.8)}
  x = x + 1
}

############################################  B - Synteny conservation and distance from promoters ############################################ 
col <- c("red", "navy", "forestgreen", "orange")
color_n = 1 # To change color between each enhancers dataset


for (enhancer in enhancers){
  sp = "mouse"
  conserv = get(paste("conserv_synteny_dist_", enhancer, sep=""))
  
  if (enhancer == "CAGE"){ # First sp
    plot(conserv[[sp]]$obs-conserv[[sp]]$simul, type="l", col=col[color_n], xaxt = "n", las=2, ylim=c(-0.1,0.1), 
         xlab="Distance from TSS (pb)", ylab="Original-Simulated \n gene-enh maintained in synteny", main="Fig.4.B")
    
  }else{lines(conserv[[sp]]$obs-conserv[[sp]]$simul, type="l", col=col[color_n])} # Add lines of other enhancers datasets
  
  # Add confidences intervals
  for (row in 1:nrow(conserv[[sp]])){
    segments(x0=row,y0=conserv[[sp]][row,"conf_low_obs"]-conserv[[sp]][row,"conf_low_simul"],
             x1=row,y1=conserv[[sp]][row,"conf_up_obs"]-conserv[[sp]][row,"conf_up_simul"],
             col=col[color_n], lwd=0.3)}
  
  color_n = color_n + 1 
}


######### Add axis and legends at the end #########
class_leg <- c("0", "500Kb", "1Mb", "1.5Mb", "2Mb", "2.5Mb")
abline(h=0, lty=2, col='black')

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.02, class_leg, xpd = TRUE)
legend("topleft", col=col, legend = enhancers, bty='n', lty=1, cex=0.8)

############################################  C - Synteny conservation and distance from promoters all species ############################################ 
species <- c("macaque", "dog", "cow", "elephant", "rabbit", "rat", "mouse", "opossum", "chicken")

col <- c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#ABD9E9","#74ADD1","#4575B4","#313695")
color_n = 1 # To change color between each enhancers dataset

conserv = get(paste("conserv_synteny_dist", enh, sep="_"))

par(lwd = 1.2) 
for (sp in species){
  if (sp == "macaque"){ # First sp
    plot(conserv[[sp]]$obs, type="l", col=col[color_n], xaxt = "n", las=2, ylim=c(0.6,1), 
         xlab="Distance from TSS (pb)", ylab="Gene-enhancer maintained in synteny ", main="Fig.4.C")
    
  }else{lines(conserv[[sp]]$obs, type="l", col=col[color_n])} # Add lines of other enhancers datasets
  
  # Add confidences intervals
  for (row in 1:nrow(conserv[[sp]])){
    segments(x0=row,y0=conserv[[sp]][row,"conf_low_obs"],x1=row,y1=conserv[[sp]][row,"conf_up_obs"], col=col[color_n], lwd=0.3)}
  
  color_n = color_n + 1 
}

par(lwd = 1) 
######### Add axis and legends at the end #########
class_leg <- c("0", "500Kb", "1Mb", "1.5Mb", "2Mb", "2.5Mb")

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.03, class_leg, xpd = TRUE)
legend("bottomleft", col=col, ncol=3, legend = species, bty='n', lty=1, cex=0.8)


dev.off()