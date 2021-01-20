############################################################################

path="/beegfs/data/necsulea/RegulatoryLandscapes/"
pathResults=paste(path, "results/co_expression_analysis/", sep="")
pathContacts="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset4/"
pathConservation="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset7/"

options(stringsAsFactors=F)
options(digits=2) ## to make files lighter

########################################################################

datasets=list()
datasets[["human"]]=c("ENCODE", "FANTOM5", "GRO-seq", "RoadmapEpigenomics")
datasets[["mouse"]]=c("FANTOM5")

synonyms=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")
names(synonyms)=c("ENCODE", "FANTOM5", "GRO-seq", "RoadmapEpigenomics")

############################################################################

for(sp in c("human", "mouse")){
  tg=setdiff(c("human", "mouse"), sp)
  
  for(dataset in datasets[[sp]]){
    print(paste(sp, dataset))

    exp=read.table(paste(pathResults, sp, "/", dataset, "/enhancer_activity_statistics.txt", sep=""), h=T, sep="\t")

    if(dataset%in%c("ENCODE", "GRO-seq", "RoadmapEpigenomics")){
      exp$id=paste(exp$chr, exp$start-1, exp$end, sep=":")
    } else{
      exp$id=paste(exp$chr, exp$start, exp$end, sep=":")
    }
    rownames(exp)=exp$id

    contacts.real=read.table(paste(pathContacts, sp, "/", synonyms[dataset], "/gene_enhancer_contacts_original_interactions.txt", sep=""), h=T)
    contacts.real=contacts.real[which(contacts.real$dist>=25e3 & contacts.real$dist<=2e6),]

    print(length(which(exp$id%in%contacts.real$enhancer)))

    exp=exp[which(exp$id%in%contacts.real$enhancer),]
    contacts.real=contacts.real[which(contacts.real$enhancer%in%exp$id),]

    ## sequence conservation

    cons=read.table(paste(pathConservation, sp, "/sequence_conservation/enhancers/",dataset, "/AlignmentStatistics_Excluding_Exons_", sp, "2", tg, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    rownames(cons)=cons[,paste("ID.",sp, sep="")]

    cons$PCUngapped=cons$FilteredUngappedLength/cons$FilteredAlignmentLength
    cons=cons[which(cons$FilteredAlignmentLength>=20),]
   
    ## median distance for enhancers

    median.dist=tapply(contacts.real$dist, as.factor(contacts.real$enhancer), median)
    mean.dist=tapply(contacts.real$dist, as.factor(contacts.real$enhancer), mean)

    exp=exp[names(median.dist),]
    cons=cons[names(median.dist),]
    cons$PCUngapped[which(is.na(cons$PCUngapped))]=0
    
    mediandist.class=cut(median.dist, breaks=seq(from=25e3, to=2e6, by=5e4), include.lowest=T)
    meandist.class=cut(mean.dist, breaks=seq(from=25e3, to=2e6, by=5e4), include.lowest=T)

    pdf(file=paste("tmp_figures/MedianContactDistance_NbSamplesExp_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$NbSamplesExp, as.factor(mediandist.class), mean), xlab="distance class", ylab="nb samples w. enhancer activity", pch=20)

    dev.off()

    pdf(file=paste("tmp_figures/MedianContactDistance_NbSamplesHighExp_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$NbSamplesHighExp, as.factor(mediandist.class), mean), xlab="distance class", ylab="nb samples w. high enhancer activity", pch=20)

    dev.off()

    pdf(file=paste("tmp_figures/MeanContactDistance_NbSamplesExp_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$NbSamplesExp, as.factor(meandist.class), mean), xlab="distance class", ylab="nb samples w. enhancer activity", pch=20)

    dev.off()

    pdf(file=paste("tmp_figures/MeanContactDistance_NbSamplesHighExp_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$NbSamplesHighExp, as.factor(meandist.class), mean), xlab="distance class", ylab="nb samples w. high enhancer activity", pch=20)

    dev.off()

     pdf(file=paste("tmp_figures/MeanContactDistance_Tau_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$Tau, as.factor(meandist.class), mean), xlab="mean distance class", ylab="mean tau", pch=20)

    dev.off()

     pdf(file=paste("tmp_figures/MeanContactDistance_TauLog_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$TauLog, as.factor(meandist.class), mean), xlab="mean distance class", ylab="mean tau log", pch=20)

    dev.off()

    
     pdf(file=paste("tmp_figures/MedianContactDistance_Tau_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$Tau, as.factor(mediandist.class), median), xlab="median distance class", ylab="median tau", pch=20)

    dev.off()

     pdf(file=paste("tmp_figures/MedianContactDistance_TauLog_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$TauLog, as.factor(mediandist.class), median), xlab="median distance class", ylab="median tau log", pch=20)

    dev.off()


    pdf(file=paste("tmp_figures/FrUngapped_",sp, "_Tau_", dataset,".pdf", sep=""), width=8, height=5)

    rho=round(cor(exp$Tau, cons$PCUngapped, method="spearman", use="complete.obs"), digits=2)
    plot(exp$Tau, cons$PCUngapped,  main=paste("rho =",rho),xlab="tau", ylab="fr ungapped", pch=20)

    dev.off()

     pdf(file=paste("tmp_figures/FrUngapped_",sp, "_NbSamplesExp_", dataset,".pdf", sep=""), width=8, height=5)

    rho=round(cor(exp$NbSamplesExp, cons$PCUngapped, method="spearman", use="complete.obs"), digits=2)
    plot(exp$NbSamplesExp, cons$PCUngapped,  xlab="tau", ylab="fr ungapped", pch=20)

    dev.off()

    
  }

}
############################################################################
