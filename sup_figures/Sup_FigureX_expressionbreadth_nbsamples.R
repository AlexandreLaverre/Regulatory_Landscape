#####################################################################

pathDatasets="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"

#####################################################################

sp="human"
assembly="hg38"

#####################################################################

expstats=read.table(paste(pathDatasets, "SupplementaryDataset6/",sp,"/ExpressionStatistics_CardosoMoreira2019.txt", sep=""), h=T, stringsAsFactors=F)

rownames(expstats)=expstats$GeneID

#####################################################################

interactions=read.table(paste(pathDatasets, "SupplementaryDataset1/",sp,"/all_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

samples=colnames(interactions)[9:dim(interactions)[2]]

interactions$nbsamples=apply(interactions[,samples], 1, function(x) length(which(!is.na(x))))

interactions$baitID=paste(interactions$chr_bait, interactions$start_bait, interactions$end_bait, sep=":")

mean.nbsamples.perbait=tapply(interactions$nbsamples, as.factor(interactions$baitID), mean)
names(mean.nbsamples.perbait)=levels(as.factor(interactions$baitID))

#####################################################################

baits=read.table(paste(pathDatasets, "SupplementaryDataset1/",sp,"/bait_coords_",assembly,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

baits$nb_genes=unlist(lapply(baits$gene_ID, function(x) length(unlist(strsplit(x, split=",")))))

baits$mean.nbsamples=mean.nbsamples.perbait[baits$ID]

baits$fac.nbsamples=cut(baits$mean.nbsamples, breaks=quantile(baits$mean.nbsamples, p=seq(from=0, to=1, length=6),  na.rm=T), include.lowest=T)

#####################################################################

baits.single.gene=baits[which(baits$nb_genes==1 & baits$gene_ID%in%rownames(expstats)),]

baits.single.gene$TauRPKM=expstats[baits.single.gene$gene_ID, "TauRPKM"]

#####################################################################

pdf("/beegfs/data/necsulea/RegulatoryLandscapesManuscript/Figures/SuppFigure7_EpressionSpecificity_NbSamplesContacts.pdf", width=6, height=5)

par(mar=c(4.5, 4.5, 1.1, 1.1))
boxplot(baits.single.gene$TauRPKM ~ baits.single.gene$fac.nbsamples, notch=T, pch=20, cex.axis=1.1)
mtext("expression specificity", line=3, side=2, cex=1.1)
mtext("mean number of samples with contacts", line=3, side=1, cex=1.1)

dev.off()

#####################################################################


