
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

##### NEW ####

path = "/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset6_regulatory_landscape_evolution/human/investigation/Gorilla/results/"
CAGE <- read.table(paste(path, "/GOrilla_results_CAGE_nb_total_process.xls", sep=""), header= T, sep="\t")
CAGE <- CAGE[which(CAGE$FDR.q.value < 5e-4),]
CAGE$regul <- grepl("regulation", CAGE$Description, fixed=T)

par(mai=c(1,3,0.1,0.5))
barplot(rev(-log10(CAGE$FDR.q.value)), horiz=T, names.arg=rev(CAGE$Description), 
        las=1, xlab="-log10(FDR)", width=0.1, cex.names=0.85,
        col=ifelse(rev(CAGE$regul) == "TRUE","mediumseagreen","grey"))

