#########################################################################
library(gsubfn)
setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")

ref_sp = "human"
 
path <- "/home/laverre/Data/Regulatory_landscape/result/"
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
if(ref_sp == "human"){rownames(expdiv)=expdiv$IDHuman}else{rownames(expdiv)=expdiv$IDMouse}

enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")

enh="CAGE"
treshold = "0.8"
data="original"

load_data <- function(enh, data){
  regland = read.table(paste(path_evol, "investigation/evolution_summary/", enh, "/", enh, "_", data, "_summary_conserv_all_", treshold, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  regland$data <- data
  common=intersect(rownames(expdiv), rownames(regland))
  
  expdiv=expdiv[common,]
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
  
  regland=regland[common,]
  regland$ratio_cons_seq = ifelse(regland$nb_total > 0, regland$nb_seq_conserv/regland$nb_total, 0)
  regland$ratio_cons_int = ifelse(regland$nb_seq_conserv > 0, regland$nb_contact_conserv/regland$nb_seq_conserv, 0)

  if (enh %in% c("ENCODE", "RoadMap")){breaks=c(1,5,10,25,50,75, max(regland$nb_total))
  breaks_names = c("1-5","6-10","11-25", "26-50", "51-75",">75")}
  
  if (enh == "CAGE"){breaks=c(1,2,5,10,15, max(regland$nb_total))
  breaks_names = c("1-2","3-5", "6-10","11-15", ">15")}
  
  if (enh == "GRO_seq"){breaks=c(1,2,5,10,15, max(regland$nb_total))
  breaks_names = c("1-2","3-5", "6-15","16-25", ">25")}
  
  regland$class_nb_contact=cut(regland$nb_total, breaks=breaks, include.lowest=T)
  regland$class_complexity = ifelse(regland$nb_total < summary(regland$nb_total)[2], "low", "high") # first quantile 
  regland$class_complexity <- factor(regland$class_complexity, levels = c("low", "high"))
  regland$class_cons_seq=cut(regland$ratio_cons_seq, breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  regland$class_cons_int=cut(regland$ratio_cons_int,  breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  
  
  return(list(expdiv, regland, breaks_names))
}


##### Regulatory landscape complexity vs Expression 
for (enh in enhancers){
  #pdf(paste(path, "Main_figures/Figure6_human_", enh, "_expression.pdf", sep=""), width=8.5, height=9.5)
  par(mai = c(0.8, 0.8, 0.2, 0.1)) # bottom, left, top, right
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  
  datas <- load_data(enh, "original")
  expdiv <- datas[[1]]
  regland <- datas[[2]]
  breaks_names <- datas[[3]]
  
  ## Expression Specificity vs. Number of enhancers
  boxplot(expdiv$TauHuman~regland$class_nb_contact, outline=F, notch=T, axes=F, xlab="", ylab="", border="navy", boxwex=0.5)
  axis(side=2)
  axis(side=1, at=1:length(breaks_names), labels=breaks_names)
  box()
  mtext("Expression specificity index", side=2, line=2.75)
  mtext("Number of contacted enhancers", side=1, line=2.75)
  
  # Expression Level vs Number of enhancers
  boxplot(log2(expdiv$Human_MeanRPKM+1)~regland$class_nb_contact, outline=F, notch=T, axes=F, xlab="", ylab="", border="navy", boxwex=0.5)
  axis(side=2)
  axis(side=1, at=1:length(breaks_names), labels=breaks_names)
  box()
  mtext("Average expression level (log2 RPKM)", side=2, line=2.75)
  mtext("Number of contacted enhancers", side=1, line=2.75)
  
  # Expression Divergence vs Number of enhancers
  boxplot(expdiv$ResidualExpressionDivergence~regland$class_nb_contact, outline=F, notch=T, axes=F, xlab="", ylab="", border="navy", boxwex=0.5)
  axis(side=2)
  axis(side=1, at=1:length(breaks_names), labels=breaks_names)
  box()
  mtext("Residual expression divergence", side=2, line=2.75)
  mtext("Number of contacted enhancers", side=1, line=2.75)
  #dev.off()

}

##### Regulatory landscape evolution 
for (enh in enhancers){
  #pdf(paste(path, "Main_figures/Figure6b_human_", enh, "_expression.pdf", sep=""), width=9, height=4)
  #par(mai = c(0.8, 0.8, 0.2, 0.1)) # bottom, left, top, right
  par(mfrow=c(1,2))
  datas <- load_data(enh, "original")
  expdiv <- datas[[1]]
  regland <- datas[[2]]
  breaks_names <- datas[[3]]
  regland_simul <- load_data(enh, "simulated")[[2]]
  selected=which(regland$nb_total>=2 & regland$nb_total<=100)

  # # Expression Divergence vs Number of conserved enhancers
  # breaks = c("0","1-25","26-40", "51-75", ">75")

  # 
  # boxplot(expdiv$ResidualExpressionDivergence[selected]~regland$class_cons_seq[selected],
  #         boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1),
  #         outline=F, notch=T, axes=F, xlab="", ylab="")
  # 
  # boxplot(expdiv$ResidualExpressionDivergence[selected]~regland$class_cons_seq[selected], add = TRUE, xaxt = "n", yaxt='n',
  #         border="darkgreen", outline=F, notch=T, boxwex=0.25, at = 1:length(breaks) - 0.2, outcex=0.2)
  # 
  # boxplot(expdiv$ResidualExpressionDivergence[selected]~regland_simul$class_cons_seq[selected], add = TRUE, xaxt = "n", yaxt='n',
  #         border="firebrick3", outline=F, notch=T, boxwex=0.25, at = 1:length(breaks) + 0.2, outcex=0.2)
  # 
  # axis(side=2)
  # axis(side=1, at=1:length(breaks), labels=breaks)
  # box()
  # mtext("Residual expression divergence", side=2, line=2.75)
  # mtext("Ratio of conserved enhancers (%)", side=1, line=2.75)
  # 
  # # Expression Divergence vs Number of conserved contacted enhancers
  # breaks = c("0","1-25","26-40", "51-75", ">75")
  # boxplot(expdiv$ResidualExpressionDivergence[selected]~regland$class_cons_int[selected],
  #         boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1),
  #         outline=F, notch=T, axes=F, xlab="", ylab="")
  # 
  # boxplot(expdiv$ResidualExpressionDivergence[selected]~regland$class_cons_int[selected], add = TRUE, xaxt = "n", yaxt='n',
  #         border="darkgreen", outline=F, notch=T, boxwex=0.25, at = 1:length(breaks) - 0.2, outcex=0.2)
  # 
  # boxplot(expdiv$ResidualExpressionDivergence[selected]~regland_simul$class_cons_int[selected], add = TRUE, xaxt = "n", yaxt='n',
  #         border="firebrick3", outline=F, notch=T, boxwex=0.25, at = 1:length(breaks) + 0.2, outcex=0.2)
  # 
  # axis(side=2)
  # axis(side=1, at=1:length(breaks), labels=breaks)
  # box()
  # mtext("Residual expression divergence", side=2, line=2.75)
  # mtext("Ratio of enhancers conserved in contact (%)", side=1, line=2.75)
  # #dev.off()
  # 
  ## Correlation between expression divergence and sequence conservation ratio
  R=cor(expdiv$ResidualExpressionDivergence[selected], regland$ratio_cons_seq[selected], method="pearson")
  rho=cor(expdiv$ResidualExpressionDivergence[selected], regland$ratio_cons_seq[selected], method="spearman")
  
  smoothScatter(expdiv$ResidualExpressionDivergence[selected], regland$ratio_cons_seq[selected], main=paste(enh, "sequence conservation"),
                xlab="ResidualExpressionDivergence", ylab="ratio_cons_seq")
  
  mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
  x=expdiv$ResidualExpressionDivergence[selected]
  y=regland$ratio_cons_seq[selected]
  abline(lm(y~x), col="red")
  
  ## Correlation between expression divergence and contacts conservation ratio
  R=cor(expdiv$ResidualExpressionDivergence[selected], regland$ratio_cons_int[selected], method="pearson")
  rho=cor(expdiv$ResidualExpressionDivergence[selected], regland$ratio_cons_int[selected], method="spearman")
  
  smoothScatter(expdiv$ResidualExpressionDivergence[selected], regland$ratio_cons_int[selected], main=paste(enh, "contact conservation"),
                xlab="ResidualExpressionDivergence", ylab="ratio_cons_int")
  
  mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
  x=expdiv$ResidualExpressionDivergence[selected]
  y=regland$ratio_cons_int[selected]
  abline(lm(y~x), col="red")
}


