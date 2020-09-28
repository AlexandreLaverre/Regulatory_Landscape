#########################################################################
setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")

ref_sp = "human"
target_sp = "mouse" 
path <- "/home/laverre/Data/Regulatory_landscape/result"
path_evol <- paste(path, "/Supplementary_dataset6_regulatory_landscape_evolution/", ref_sp, "/", sep="")

compute.tau <- function(exp){
  if(max(exp)==0){
    return(NA)
  }
  
  n=length(exp)
  newexp=exp/max(exp)

  tau=sum(1-newexp)/(n-1)

  return(tau)
}

#############################################################################
expdiv=read.table("ExpressionDivergence_CardosoMoreira2019.txt", h=T, stringsAsFactors=F, sep="\t")

#expdiv=read.table("../../results/expression_divergence/ExpressionDivergence_CardosoMoreira2019.txt", h=T, stringsAsFactors=F, sep="\t")
rownames(expdiv)=expdiv$IDHuman

#regland=read.table("human_merged_to_mouse_summary_by_gene.txt", h=T, stringsAsFactors=F, sep="\t")
#regland=read.table("../../results/landscape_conservation/human2mouse_conservation_by_gene_0.4.txt", h=T, stringsAsFactors=F, sep="\t")
enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")


enh="ENCODE"
treshold = "0.5"
file = paste("/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset6_regulatory_landscape_evolution/human/", enh, "_original_evolution_summary_all_100000.txt", sep="")
regland = read.table(file, h=T, stringsAsFactors=F, sep="\t", row.names = 1)

#regland = read.table(paste(path_evol, enh, "_original_summary_conserv_", treshold, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
#regland_simul = read.table(paste(path_evol, enh, "_simulated_summary_conserv_", treshold, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)

common=intersect(rownames(expdiv), rownames(regland))
common=intersect(common, rownames(regland_simul))

expdiv=expdiv[common,]
regland=regland[common,]
regland_simul=regland_simul[common,]

regland$ratio_cons_int=regland$nb_contact_conserv/regland$nb_seq_conserv
regland$ratio_cons_seq=regland$nb_seq_conserv/regland$nb_total

if (enh == "ENCODE"){breaks=c(1,5,10,25,50,75, max(regland$nb_total))
breaks_names = c("1-5","6-10","11-25", "26-50", "51-75",">75")

}else if (enh == "CAGE"){breaks=c(1,2,5,10,15, max(regland$nb_total))
breaks_names = c("1-2","3-5", "6-10","11-15", ">15")}

regland$class_nb_contact=cut(regland$nb_total, breaks=breaks, include.lowest=T)

regland_simul$ratio_cons_int=regland_simul$nb_contact_conserv/regland_simul$nb_seq_conserv
regland_simul$ratio_cons_seq=regland_simul$nb_seq_conserv/regland_simul$nb_total
regland_simul$class_nb_contact=cut(regland_simul$nb_total, breaks=breaks, include.lowest=T)


#########################################################################

minRPKM=1

samples.human=setdiff(grep("^Human_", colnames(expdiv), value=T), "Human_MeanRPKM")
samples.mouse=setdiff(grep("^Mouse_", colnames(expdiv), value=T), "Mouse_MeanRPKM")

tissues.human=unlist(lapply(samples.human, function(x) unlist(strsplit(x, split="_"))[2]))
tissues.mouse=unlist(lapply(samples.mouse, function(x) unlist(strsplit(x, split="_"))[2]))

stage.human=unlist(lapply(samples.human, function(x) unlist(strsplit(x, split="_"))[3]))
stage.mouse=unlist(lapply(samples.mouse, function(x) unlist(strsplit(x, split="_"))[3]))

expdiv$NbSamplesMouse=apply(expdiv[,samples.mouse], 1, function(x) length(which(x>=minRPKM)))
expdiv$NbSamplesHuman=apply(expdiv[,samples.human], 1, function(x) length(which(x>=minRPKM)))

expdiv$TauHuman=apply(expdiv[,samples.human],1, function(x) compute.tau(as.numeric(x)))
expdiv$TauMouse=apply(expdiv[,samples.mouse],1, function(x) compute.tau(as.numeric(x)))

col.tissues=c("navy", "steelblue", "indianred", "seagreen", "orange")
names(col.tissues)=c("Brain", "Cerebellum", "Heart", "Kidney", "Liver")

############################ RPKM SHH #############################################

shh.human=as.numeric(expdiv["ENSG00000164690", samples.human])
shh.mouse=as.numeric(expdiv["ENSG00000164690", samples.mouse])


#pdf(file="figures/RPKM_Shh.pdf", width=8.5, height=7)
par(mfrow=c(2,1))
par(mar=c(3.1,4.1,1.1,1.1))
b=barplot(shh.human, col=col.tissues[tissues.human], border=col.tissues[tissues.human], axes=F)
axis(side=2, cex=0.95)
legend("topright", legend=c("forebrain", "cerebellum", "kidney", "liver", "heart"), fill=col.tissues[c("Brain", "Cerebellum",  "Kidney", "Liver","Heart")], border=col.tissues[c("Brain", "Cerebellum",  "Kidney", "Liver","Heart")], inset=0.01, bty="n")
mtext("RPKM, human samples", side=2, line=3)
mtext(stage.human, side=1, at=b, las=2, line=0.5, cex=0.8)

b=barplot(shh.mouse, col=col.tissues[tissues.mouse], border=col.tissues[tissues.mouse], axes=F)
axis(side=2, cex=0.95)
mtext("RPKM, mouse samples", side=2, line=3)
mtext(stage.mouse, side=1, at=b, las=2, line=0.5, cex=0.8)
#dev.off()

############################ Relative expression profile SHH ###########################################

rel.shh.human=shh.human/sum(shh.human)
rel.shh.mouse=shh.mouse/sum(shh.mouse)

#pdf(file="figures/RelativeExpressionLevel_Shh.pdf", width=8.5, height=7)
par(mfrow=c(2,1))
par(mar=c(3.1,4.1,1.1,1.1))
b=barplot(rel.shh.human, col=col.tissues[tissues.human], border=col.tissues[tissues.human], axes=F)
axis(side=2, cex=0.95)
legend("topright", legend=c("forebrain", "cerebellum", "kidney", "liver", "heart"), fill=col.tissues[c("Brain", "Cerebellum",  "Kidney", "Liver","Heart")], border=col.tissues[c("Brain", "Cerebellum",  "Kidney", "Liver","Heart")], inset=0.01, bty="n")
mtext("relative expression level, human", side=2, line=3)
mtext(stage.human, side=1, at=b, las=2, line=0.5, cex=0.8)

b=barplot(rel.shh.mouse, col=col.tissues[tissues.mouse], border=col.tissues[tissues.mouse], axes=F)
axis(side=2, cex=0.95)
mtext("relative expression level, mouse", side=2, line=3)
mtext(stage.mouse, side=1, at=b, las=2, line=0.5, cex=0.8)
#dev.off()

#########################################################################

## expression divergence vs. mean expression level

#pdf("figures/ExpressionDivergence_MeanExpressionLevel.pdf", width=7, height=7)

smoothScatter(log2(expdiv$MeanRPKM+1), expdiv$ExpressionDivergence, xlab="Mean expression level (log2-transformed RPKM)", ylab="Expression divergence")
rho=cor(log2(expdiv$MeanRPKM+1), expdiv$ExpressionDivergence, method="spearman")
R=cor(log2(expdiv$MeanRPKM+1), expdiv$ExpressionDivergence, method="pearson")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
abline(lm(expdiv$ExpressionDivergence~log2(expdiv$MeanRPKM+1)), col="red")

#dev.off()
       
#########################################################################

## expression specificity vs. regulatory landscape

#pdf(file="figures/ExpressionSpecificity_vs_NbContacts.pdf", width=8, height=6)
boxplot(expdiv$TauHuman~regland$class_nb_contact, outline=F, notch=T, axes=F, xlab="", ylab="")
axis(side=2)
axis(side=1, at=1:length(breaks_names), labels=breaks_names)
box()
mtext("Expression specificity index", side=2, line=2.75)
mtext("Number of contacted regions", side=1, line=2.75)
#dev.off()

#########################################################################

#pdf(file="figures/MeanExpressionLevel_vs_NbContacts.pdf", width=8, height=6)
boxplot(log2(expdiv$Human_MeanRPKM+1)~regland$class_nb_contact, outline=F, notch=T, axes=F, xlab="", ylab="")
axis(side=2)
axis(side=1, at=1:length(breaks_names), labels=breaks_names)
box()
mtext("Average expression level (log2 RPKM)", side=2, line=2.75)
mtext("Number of contacted regions", side=1, line=2.75)
#dev.off()

#########################################################################

#### expression divergence as a function of the number of contacts

#pdf("figures/ExpressionDivergence_NbContacts_2classes.pdf",width=4, height=5)
boxplot(expdiv$ExpressionDivergence~regland$class_nb_contact, outline=F, notch=T,
        xlab="Total number of contacts", ylab="Expression divergence", names=breaks_names, boxwex = 0.7)
#dev.off()


#pdf("figures/ResidualExpressionDivergence_NbContacts_Class.pdf")
boxplot(expdiv$ResidualExpressionDivergence~regland$class_nb_contact, outline=F, notch=T, boxwex = 0.7,
        xlab="Total number of contacts", ylab="Residual expression divergence", names=breaks_names)
#dev.off()

###########################################################################
library(RColorBrewer)
if (enh="ENCODE"){min}
selected=which(regland$nb_total>=10 & regland$nb_total<=100)

#pdf("figures/ExpressionDivergence_RatioConservedContacts.pdf")
R=cor(regland$ratio_cons_int[selected],expdiv$ExpressionDivergence[selected], method="pearson")
rho=cor(regland$ratio_cons_int[selected],expdiv$ExpressionDivergence[selected], method="spearman")

smoothScatter(100*regland$ratio_cons_int[selected], expdiv$ExpressionDivergence[selected], colramp = colorRampPalette(brewer.pal(8, "Greens")),
              xlab="% conserved contacts", ylab="Expression divergence")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
x=(100*regland$ratio_cons_int[selected])
y=expdiv$ExpressionDivergence[selected]
abline(lm(y~x), col="black")

# Simul
R=cor(regland_simul$ratio_cons_int[selected],expdiv$ExpressionDivergence[selected], method="pearson")
rho=cor(regland_simul$ratio_cons_int[selected],expdiv$ExpressionDivergence[selected], method="spearman")

smoothScatter(100*regland_simul$ratio_cons_int[selected], expdiv$ExpressionDivergence[selected], 
              colramp = colorRampPalette(c(rgb(1, 1, 1, 0), rgb(1, 0, 0, 1)), alpha = TRUE),
              xlab="% conserved contacts", ylab="Expression divergence", xlim=c(0,100), add=TRUE)

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
x=(100*regland_simul$ratio_cons_int[selected])
y=expdiv$ExpressionDivergence[selected]
abline(lm(y~x), col="red")
#dev.off()


#pdf("figures/ExpressionDivergence_RatioConservedSequences.pdf")
R=cor(regland$ratio_cons_seq[selected],expdiv$ExpressionDivergence[selected], method="pearson")
rho=cor(regland$ratio_cons_seq[selected],expdiv$ExpressionDivergence[selected], method="spearman")

smoothScatter(100*regland$ratio_cons_seq[selected], expdiv$ExpressionDivergence[selected],  xlab="% conserved sequences", ylab="Expression divergence")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
y=expdiv$ExpressionDivergence[selected]
x=(100*regland$ratio_cons_seq[selected])
abline(lm(y~x), col="red")
#dev.off()

###########################################################################

#pdf("figures/ResidualExpressionDivergence_RatioConservedContacts.pdf")
R=cor(regland$ratio_cons_int[selected],expdiv$ResidualExpressionDivergence[selected], method="pearson")
rho=cor(regland$ratio_cons_int[selected],expdiv$ResidualExpressionDivergence[selected], method="spearman")

smoothScatter(100*regland$ratio_cons_int[selected], expdiv$ResidualExpressionDivergence[selected],  xlab="% conserved contacts", ylab="Residual expression divergence")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
x=(100*regland$ratio_cons_int[selected])
y=expdiv$ResidualExpressionDivergence[selected]
abline(lm(y~x), col="red")
#dev.off()


#pdf("figures/ResidualExpressionDivergence_RatioConservedSequences.pdf")
R=cor(regland$ratio_cons_seq[selected],expdiv$ResidualExpressionDivergence[selected], method="pearson")
rho=cor(regland$ratio_cons_seq[selected],expdiv$ResidualExpressionDivergence[selected], method="spearman")

smoothScatter(100*regland$ratio_cons_seq[selected], expdiv$ResidualExpressionDivergence[selected],  xlab="% conserved sequences", ylab="Residual expression divergence")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
y=expdiv$ResidualExpressionDivergence[selected]
x=(100*regland$ratio_cons_seq[selected])
abline(lm(y~x), col="red")
#dev.off()


###########################################################################

