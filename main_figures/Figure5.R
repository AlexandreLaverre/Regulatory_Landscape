#########################################################################################################################
source("parameters.R") ## pathFiguress are defined based on the user name

library(ape)
library(vioplot)

# Choose species within : human ; mouse
ref_sp = "human" 
# Choose genes within : all ; dvpt ; other
selected_genes = "all"

target_sp = "mouse"

enhancers <- c("FANTOM5", "ENCODE")
if (ref_sp == "human"){enh_names <- c(enhancers, "RoadMap\nEpigenomics", "GRO-seq"); enhancers <- c(enhancers, "RoadmapEpigenomics", "FOCS_GRO_seq")}

load(paste(pathFigures, "Fig5_", ref_sp, "_", selected_genes, "_genes.Rdata", sep=""))

col <- c("red", "navy", "forestgreen", "orange")
#########################################################################################################################
if(ref_sp == "human"){pdf_name="Figure5_bis.pdf"}else{pdf_name="Sup_Figure14.pdf"}

pdf(paste(pathFigures, pdf_name, sep=""), width=8.5, height=8)
par(mai = c(0.8, 0.8, 0.5, 0.1)) # bottom, left, top, right
layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))

############################################   A - Global contact conservation ############################################ 
if (ref_sp == "human"){YMAX=25}else{YMAX=30}

par(lwd = 1.5) 
bar <- barplot(conserv_global$data, border=rep(c("darkgreen", "firebrick3", "black"),4), beside=T, space = c(0, 0.1, 0),
               ylim=c(0,YMAX), col="white", main="", ylab="Conserved contact (%)", las=2)

arrows(x0=bar,y0=conserv_global$conf_up,y1=conserv_global$conf_low,angle=90,code=3,length=0.05)
text(conserv_global$n_total, x=bar, y=conserv_global$conf_up+2, cex=0.7)

axis(side=1, at=c(1,4.2,7.4,10.4), labels=enh_names, mgp=c(3, 1.2, 0), cex.axis=0.8)
#text(c(1,4.2,7.4,10.4), par("usr")[3]-0.005, labels = enh_names, pos = 1, xpd = TRUE)  
legend("topright", legend = c("Original", "Simulated"), border=c("darkgreen", "firebrick3"), fill="white", bty='n')
mtext("A", side=3, line=1, at=0.1, font=2, cex=1.2)

############################################  B - Contact conservation by distance from TSS ############################################ 
if (ref_sp == "human"){YLIM=c(-3,8)}else{YLIM=c(-5, 20)}
color_n = 1 # To change color between each enhancers dataset
class_leg <- c("0",  "0.5",  "1", "1.5", "2")

par(lwd = 0.7) 
for (enh in enhancers){
  conserv = get(paste("conserv_dist", enh, sep="_"))
  if (enh == "FANTOM5"){ # First enhancers dataset
    plot(log((conserv$result-conserv$simul)/conserv$simul), type="l", col=col[color_n], xaxt = "n", ylim=YLIM,
         xlab="Distance from TSS (Mb)", ylab="Excess of contact conservation\n from expected (log)", main="", las=2)
  }else{lines(log((conserv$result-conserv$simul)/conserv$simul), type="l", col=col[color_n])} # Add lines of other enhancers datasets
  
  for (row in 1:nrow(conserv)){
    segments(x0=row,y0=log((conserv[row,]$obs_conflow-conserv[row,]$simul_conflow)/conserv[row,]$simul_conflow),
             x1=row,y1=log((conserv[row,]$obs_confup-conserv[row,]$simul_confup)/conserv[row,]$simul_confup),
             col=col[color_n], lwd=0.3)}
  
  color_n = color_n + 1
  
}


axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
legend("topleft", col=col, legend = enh_names, bty='n', lty=1, cex=0.8)
mtext("B", side=3, line=1, at=0.1, font=2, cex=1.2)

############################################  C - Contact conservation by nb samples ############################################ 
color_n = 1
par(lwd = 1.5)
enh = "ENCODE"
conserv = get(paste("conserv_sample", enh, sep="_"))

bar <- barplot(result ~ data+class, beside=T, data=conserv, space = c(0.1, 0.6), ylim=c(0,60),
               border=c("darkgreen", "firebrick1"), col="white", ylab="Conserved contact (%)", xlab="Number of sample", 
               main="")

arrows(x0=bar,y0=conserv[c(1,6,2,7,3,8,4,9,5,10),]$confup,y1=conserv[c(1,6,2,7,3,8,4,9,5,10),]$conflow,angle=90,code=3,length=0.05)
legend("topright", fill="white", border=c("darkgreen", "firebrick"), legend = c("Original", "Simulated"), bty='n', cex=0.9)

text(conserv[c(1,6,2,7,3,8,4,9,5,10),]$count,x = bar, y=conserv[c(1,6,2,7,3,8,4,9,5,10),]$confup+3, cex=0.7)
mtext("C", side=3, line=1, at=0.1, font=2, cex=1.2)

############################################ D - Contact conservation by similar samples ####################################
if (ref_sp == "human"){YMAX=30; side="topright"}else{YMAX=17; side="topleft"}

par(lwd = 1.5) 
conserv = get(paste("conserv_similar_sample", enh, sep="_"))

bar <- barplot(conserv$data, border=rep(c("darkgreen", "firebrick3", "black"),4), beside=T, space = c(0, 0.1, 0), las=2,
               ylim=c(0,YMAX), col="white", main="", ylab="Conserved contact (%)")
arrows(x0=bar,y0=conserv$conf_up,y1=conserv$conf_low,angle=90,code=3,length=0.05)
text(conserv$n_total, x=bar, y=conserv$conf_up+2, cex=0.8)

cells = c("Pre-adipocytes", "ESC", "Bcell")
text(c(1,4.2,7.4), par("usr")[3]-0.5, labels = cells, pos = 1, xpd = TRUE)  
legend(side, legend = c("Original", "Simulated"), border=c("darkgreen", "firebrick3"), fill="white", bty='n')
mtext("D", side=3, line=1, at=0.1, font=2, cex=1.2)

dev.off()
