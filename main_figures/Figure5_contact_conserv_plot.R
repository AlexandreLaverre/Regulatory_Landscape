############################################################ PLOT FIGURE 5 ###############################################
path <- "/home/laverre/Data/Regulatory_landscape/result/Main_figures/"
load(paste(path, "Fig5_human_corrected.Rdata", sep=""))


pdf(paste(path, "Figure5_human_corrected.pdf", sep=""), width=8.5, height=8)
par(mai = c(0.8, 0.8, 0.5, 0.1)) # bottom, left, top, right
layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))

enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")
    
############################################   A - Global contact conservation ############################################ 
par(lwd = 1.5) 
bar <- barplot(conserv_global$data, border=rep(c("darkgreen", "firebrick3", "black"),4), beside=T, space = c(0, 0.1, 0),
               ylim=c(0,26), col="white", main="Fig 5.A.", ylab="Conserved contact (%)", las=2)
arrows(x0=bar,y0=conserv_global$conf_up,y1=conserv_global$conf_low,angle=90,code=3,length=0.05)
text(conserv_global$n_total, x=bar, y=conserv_global$conf_up+2, cex=0.7)

text(c(1,4.2,7.4,10.4), par("usr")[3]-0.005, labels = enhancers, pos = 1, xpd = TRUE)  
legend("topright", legend = c("Original", "Simulated"), border=c("darkgreen", "firebrick3"), fill="white", bty='n')


############################################  B - Contact conservation by distance from TSS ############################################ 
col <- c("red", "navy", "forestgreen", "orange")
color_n = 1 # To change color between each enhancers dataset
class_leg <- c("0",  "500kb",  "1Mb", "1.5Mb", "2Mb", "2.5Mb", "3Mb")

par(lwd = 0.7) 
for (enh in enhancers){
  conserv = get(paste("conserv_dist", enh, sep="_"))
  if (enh == "CAGE"){ # First enhancers dataset
    plot(conserv$result-conserv$simul, type="l", col=col[color_n], xaxt = "n", ylim=c(0,17),
         xlab="Distance from TSS (pb)", ylab="Original - Simulated \n conserved contact (%)", main="Fig 5.B.", las=2)
  }else{lines(conserv$result-conserv$simul, type="l", col=col[color_n])} # Add lines of other enhancers datasets
  
  for (row in 1:nrow(conserv)){
    segments(x0=row,y0=conserv[row,]$obs_conflow-conserv[row,]$simul_conflow,x1=row,y1=conserv[row,]$obs_confup-conserv[row,]$simul_confup, col=col[color_n], lwd=0.3)}
  
  color_n = color_n + 1

}

axis(1, at=seq(1,61,10), labels=F)
text(seq(1,61,10),par("usr")[3]-1, class_leg, xpd = TRUE)
legend("topleft", col=col, legend = enhancers, bty='n', lty=1, cex=0.8)


############################################  C - Contact conservation by nb samples ############################################ 
color_n = 1
par(lwd = 1.5)
enh = "ENCODE"
conserv = get(paste("conserv_sample", enh, sep="_"))

bar <- barplot(result ~ data+class, beside=T, data=conserv, space = c(0.1, 0.6), ylim=c(0,60),
               border=c("darkgreen", "firebrick1"), col="white", ylab="Conserved contact (%)", xlab="Number of sample", 
               main="Fig 5.C.")

arrows(x0=bar,y0=conserv[c(1,6,2,7,3,8,4,9,5,10),]$confup,y1=conserv[c(1,6,2,7,3,8,4,9,5,10),]$conflow,angle=90,code=3,length=0.05)
legend("topright", fill="white", border=c("darkgreen", "firebrick"), legend = c("Original", "Simulated"), bty='n', cex=0.9)

text(conserv[c(1,6,2,7,3,8,4,9,5,10),]$count,x = bar, y=conserv[c(1,6,2,7,3,8,4,9,5,10),]$confup+3, cex=0.7)

############################################ D - Contact conservation by similar samples ####################################
par(lwd = 1.5) 
conserv = get(paste("conserv_similar_sample", enh, sep="_"))

bar <- barplot(conserv$data, border=rep(c("darkgreen", "firebrick3", "black"),4), beside=T, space = c(0, 0.1, 0), las=2,
               ylim=c(0,26), col="white", main="Fig 5.D.", ylab="Conserved contact (%)")
arrows(x0=bar,y0=conserv$conf_up,y1=conserv$conf_low,angle=90,code=3,length=0.05)
text(conserv$n_total, x=bar, y=conserv$conf_up+2, cex=0.8)

cells = c("Pre-adipocytes", "ESC", "Bcell")
text(c(1,4.2,7.4), par("usr")[3]-0.5, labels = cells, pos = 1, xpd = TRUE)  
legend("topright", legend = c("Original", "Simulated"), border=c("darkgreen", "firebrick3"), fill="white", bty='n')

dev.off()