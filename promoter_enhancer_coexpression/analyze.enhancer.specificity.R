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

    cons=read.table(paste(pathConservation, sp, "/sequence_conservation/enhancers/",synonyms[dataset], "/AlignmentStatistics_Excluding_Exons_", sp, "2", tg, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    rownames(cons)=cons[,paste("ID.",sp, sep="")]

    cons$PCUngapped=cons$FilteredUngappedLength/cons$FilteredAlignmentLength
    cons=cons[which(cons$FilteredAlignmentLength>=20),]
   
    ## median distance for enhancers

    median.dist=tapply(contacts.real$dist, as.factor(contacts.real$enhancer), median)
    mean.dist=tapply(contacts.real$dist, as.factor(contacts.real$enhancer), mean)

    mediandist.class=cut(median.dist, breaks=seq(from=25e3, to=2e6, by=5e4), include.lowest=T)
    meandist.class=cut(mean.dist, breaks=seq(from=25e3, to=2e6, by=5e4), include.lowest=T)

    exp=exp[names(median.dist),]

    ## pc ungapped
    pcungap=rep(NA, dim(exp)[1])
    names(pcungap)=rownames(exp)

    pcungap[intersect(rownames(exp), rownames(cons))]=cons[intersect(rownames(exp), rownames(cons)), "PCUngapped"]

    tauclass=cut(exp$Tau, breaks=quantile(exp$Tau, p=seq(from=0, to=1, length=11), na.rm=T), include.lowest=T)

    ## tau vs distance & sequence conservation
    
    ## pc ungapped as a function of tau

    pdf(file=paste("tmp_figures/TauClass_PCUngapped_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(pcungap, tauclass, mean, na.rm=T), xlab="Tau class", ylab="mean fraction ungapped sequence", pch=20)

    dev.off()
  
   
    pdf(file=paste("tmp_figures/MedianContactDistance_Tau_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$Tau, as.factor(mediandist.class), mean), xlab="median distance class", ylab="mean tau", pch=20)

    dev.off()

    ## nb samples with high exp

    nbsamplesclass=cut(exp$NbSamplesTPM1, breaks=quantile(exp$NbSamplesTPM1, p=seq(from=0, to=1, length=11), na.rm=T), include.lowest=T)

     pdf(file=paste("tmp_figures/NbSamplesExp_PCUngapped_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(pcungap, nbsamplesclass, mean, na.rm=T), xlab="nbsamples class", ylab="mean fraction ungapped sequence", pch=20)

    dev.off()
  
   
    pdf(file=paste("tmp_figures/MedianContactDistance_NbSamplesExpTPM1_",sp, "_", dataset,".pdf", sep=""), width=8, height=5)

    plot(tapply(exp$NbSamplesTPM1, as.factor(mediandist.class), mean), xlab="median distance class", ylab="mean nb samples TPM 1", pch=20)

    dev.off()
    
  }

}
############################################################################
