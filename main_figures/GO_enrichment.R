setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/GOrilla/")

sp = "human"
#CAGE <- read.table(paste(sp, "/GO_", sp,"_gene_syntenie_sup85.xls", sep=""), header= T, sep="\t")
CAGE <- read.table(paste(sp, "/CAGE_sup15.xls", sep=""), header= T, sep="\t")

CAGE <- CAGE[which(CAGE$FDR.q.value < 2.5e-6 & CAGE$P.value < 2.5e-9),]
#CAGE$regul <- grepl("regulation", CAGE$Description, fixed=T)
7.55e-13
5.61e-10

par(mai=c(1,3,0.1,0.5))
barplot(rev(-log10(CAGE$FDR.q.value)), horiz=T, names.arg=rev(CAGE$Description), 
        las=1, xlab="-log10(FDR)", width=0.1, cex.names=0.85)


#,col=ifelse(rev(CAGE$regul) == "TRUE","mediumseagreen","grey"))
#legend("bottomright", legend=c("regulation", "other"), fill=c("mediumseagreen", "grey"), bty="n")

