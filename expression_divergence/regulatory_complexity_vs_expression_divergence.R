#########################################################################
library(gsubfn)
setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")
path <- "/home/laverre/Data/Regulatory_landscape/result/"

enh="CAGE"        # CAGE or ENCODE
treshold = "0.5"  # alignment score
data="simulated"   # original or simulated
sample="all"      # all ; pre_adipo ; Bcell or ESC

############################################################## Data ############################################################## 

expdiv=read.table("ExpressionDivergence_CardosoMoreira2019.txt", h=T, stringsAsFactors=F, sep="\t")

compute.tau <- function(exp){
  if(max(exp)==0){
    return(NA)
  }
  
  n=length(exp)
  newexp=exp/max(exp)
  
  tau=sum(1-newexp)/(n-1)
  
  return(tau)
}

samples.human=setdiff(grep("^Human_", colnames(expdiv), value=T), "Human_MeanRPKM")
samples.mouse=setdiff(grep("^Mouse_", colnames(expdiv), value=T), "Mouse_MeanRPKM")

expdiv$TauHuman=apply(expdiv[,samples.human],1, function(x) compute.tau(as.numeric(x)))
expdiv$TauMouse=apply(expdiv[,samples.mouse],1, function(x) compute.tau(as.numeric(x)))


### Mouse score
ref_sp = "mouse"
path_evol <- paste(path, "/Supplementary_dataset6_regulatory_landscape_evolution/", ref_sp, "/", sep="")

regland = read.table(paste(path_evol, "investigation/evolution_summary/", enh, "_", data, "_summary_conserv_", sample, "_", treshold, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
rownames(expdiv)=expdiv$IDMouse
common=intersect(rownames(expdiv), rownames(regland))

expdiv=expdiv[common,]
regland=regland[common,]
expdiv$mouse_zscore <- (regland$nb_total-mean(regland$nb_total)) / sd(regland$nb_total)
expdiv <- expdiv[which(expdiv$mouse_zscore < 10),]

### Human score
ref_sp = "human"
path_evol <- paste(path, "/Supplementary_dataset6_regulatory_landscape_evolution/", ref_sp, "/", sep="")
regland = read.table(paste(path_evol, "investigation/evolution_summary/", enh, "_", data, "_summary_conserv_", sample, "_", treshold, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
rownames(expdiv)=expdiv$IDHuman
common=intersect(rownames(expdiv), rownames(regland))

expdiv=expdiv[common,]
regland=regland[common,]
expdiv$human_zscore <- (regland$nb_total-mean(regland$nb_total)) / sd(regland$nb_total)
expdiv <- expdiv[which(expdiv$human_zscore < 10),]

expdiv$delta_complex <- abs(expdiv$mouse_zscore-expdiv$human_zscore)
############################################################## Results ############################################################## 
#pdf("/home/laverre/Data/Regulatory_landscape/result/evol_complexity_all.pdf")

## Correlation between regulatory landscapes complexity
R=cor(expdiv$human_zscore,expdiv$mouse_zscore, method="pearson")
rho=cor(expdiv$human_zscore,expdiv$mouse_zscore, method="spearman")

smoothScatter(expdiv$human_zscore, expdiv$mouse_zscore, main=paste(enh, "in", sample, "samples"),
              xlab="Human complexity zscore", ylab="Mouse complexity zscore")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
x=expdiv$human_zscore
y=expdiv$mouse_zscore
abline(lm(y~x), col="red")

## Correlation between RPKM and complexity
R=cor(expdiv$delta_complex,log2(expdiv$MeanRPKM+1), method="pearson")
rho=cor(expdiv$delta_complex,log2(expdiv$MeanRPKM+1), method="spearman")

smoothScatter(expdiv$delta_complex, log2(expdiv$MeanRPKM+1), main=paste(enh, "in", sample, "samples"),
              xlab="Delta complexity", ylab="MeanRPKM")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
x=expdiv$delta_complex
y=log2(expdiv$MeanRPKM+1)
abline(lm(y~x), col="red")

## Correlation between Specificity and complexity
R=cor(expdiv$delta_complex,expdiv$mouse_zscore, method="pearson")
rho=cor(expdiv$delta_complex,expdiv$mouse_zscore, method="spearman")

smoothScatter(expdiv$delta_complex, expdiv$mouse_zscore, main=paste(enh, "in", sample, "samples"),
              xlab="Delta complexity", ylab="mouse_zscore")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
x=expdiv$delta_complex
y=expdiv$mouse_zscore
abline(lm(y~x), col="red")


## Correlation between Divergence Expression and divergence of complexity
R=cor(expdiv$delta_complex,expdiv$ExpressionDivergence, method="pearson")
rho=cor(expdiv$delta_complex,expdiv$ExpressionDivergence, method="spearman")

smoothScatter(expdiv$delta_complex, expdiv$ExpressionDivergence, main=paste(enh, "in", sample, "samples"),
              xlab="Delta complexity", ylab="Expression Divergence")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
x=expdiv$delta_complex
y=expdiv$ExpressionDivergence
abline(lm(y~x), col="red")

## Correlation between Residual Divergence Expression and divergence of complexity 
R=cor(expdiv$delta_complex,expdiv$ResidualExpressionDivergence, method="pearson")
rho=cor(expdiv$delta_complex,expdiv$ResidualExpressionDivergence, method="spearman")

smoothScatter(expdiv$delta_complex, expdiv$ResidualExpressionDivergence, main=paste(enh, "in", sample, "samples"),
              xlab="Delta complexity", ylab="Residual Expression Divergence")

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
x=expdiv$delta_complex
y=expdiv$ResidualExpressionDivergence
abline(lm(y~x), col="red")

dev.off()