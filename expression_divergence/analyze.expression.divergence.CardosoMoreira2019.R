#########################################################################

expdiv=read.table("../../results/expression_divergence/ExpressionDivergence_CardosoMoreira2019.txt", h=T, stringsAsFactors=F, sep="\t")
rownames(expdiv)=expdiv$IDHuman

expdiv$MaxExp=apply(as.matrix(expdiv[,-c(1,2,3)]), 1, max)
lm1=lm(expdiv$ExpressionDivergence~expdiv$MaxExp)
expdiv$ResidualExpressionDivergence

regland=read.table("../../results/landscape_conservation/human2mouse_conservation_by_gene_0.4.txt", h=T, stringsAsFactors=F, sep="\t")
rownames(regland)=regland$TSS

common=intersect(rownames(expdiv), rownames(regland))
expdiv=expdiv[common,]
regland=regland[common,]

regland$ratio_cons_int=regland$nb_int_conserv/regland$nb_contact
regland$ratio_cons_seq=regland$nb_seq_conserv/regland$nb_contact

#########################################################################

## expression divergence as a function of the number of contacts

fac.nbcontacts=cut(regland$nb_contact, breaks=c(1,10, 20, 30, 40, 50, max(regland$nb_contact)), include.lowest=T)

pdf("figures/ExpDivergence_NbContacts_Class.pdf")
boxplot(expdiv$ExpressionDivergence~fac.nbcontacts, outline=F, notch=T, xlab="Total number of contacts", ylab="Expression divergence")
dev.off()


fac2.nbcontacts=cut(regland$nb_contact, breaks=c(1,10, max(regland$nb_contact)), include.lowest=T)

pdf("figures/ExpDivergence_NbContacts_2classes.pdf", width=4.5, height=6)
boxplot(expdiv$ExpressionDivergence~fac2.nbcontacts, outline=F, notch=T, xlab="Total number of contacts", ylab="Expression divergence", names=c("1 to 10", ">10"), boxwex=0.5, cex.axis=1.25, cex.lab=1.25)
dev.off()


pdf("figures/ExpDivergence_NbContacts.pdf")
smoothScatter(regland$nb_contact, expdiv$ExpressionDivergence,  xlab="Total number of contacts", ylab="Expression divergence")
dev.off()


pdf("figures/ExpDivergence_RatioConservedContacts.pdf")
smoothScatter(100*regland$ratio_cons_int, expdiv$ExpressionDivergence,  xlab="% conserved contacts", ylab="Expression divergence")
dev.off()

#########################################################################

#########################################################################




