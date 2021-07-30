###############################################################################

pathFinalData="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
pathContactedEnhancers=paste(pathFinalData, "SupplementaryDataset4/", sep="")
pathNeighborEnhancers="/beegfs/data/necsulea/RegulatoryLandscapes/results/neighbor_enhancers/"

release=94

###############################################################################

enhancer.datasets=list()
enhancer.datasets[["human"]]=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")
enhancer.datasets[["mouse"]]=c("ENCODE", "FANTOM5")

###############################################################################

for(sp in c("human", "mouse")){
  for(enh in enhancer.datasets[[sp]]){

    print(paste(sp, enh))

    n=read.table(paste(pathNeighborEnhancers, sp, "/", enh, "/neighbor_enhancers_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    g=grep("^chr", n$EnhancerID)

    if(length(g)==0){
      n$EnhancerID=paste("chr", n$EnhancerID, sep="")
    }
    
    c=read.table(paste(pathContactedEnhancers, sp, "/", enh, "/gene_enhancer_contacts_original_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    c=c[which(c$min_dist_baitedTSS>=25000 & c$min_dist_baitedTSS<=2e6),]

    ## select common genes
    common.genes=intersect(n$GeneID, c$gene)

    n=n[which(n$GeneID%in%common.genes),]
    c=c[which(c$gene%in%common.genes),]

    ## nb enhancers per gene

    n.nbenh=tapply(n$EnhancerID, as.factor(n$GeneID), function(x) length(unique(x)))
    c.nbenh=tapply(c$enhancer, as.factor(c$gene), function(x) length(unique(x)))

    pdf(file=paste("figures/Comparison_NbEnhancers_",sp,"_",enh,".pdf",sep=""), width=6, height=6)
    par(mar=c(4.5, 4.5, 2.1, 1.1))

    plot(c.nbenh, n.nbenh, pch=20, xlab="nb contacted enhancers", ylab="nb neighboring enhancers",main="")

    R=round(cor(n.nbenh, c.nbenh), digits=2)
    rho=round(cor(n.nbenh, c.nbenh, method="spearman"), digits=2)

    mtext(paste("R = ", R, " rho = ", rho, sep=""), line=0.5, side=3)

    abline(0, 1, col="red")
    dev.off()

    ## compare distance distributions

    pdf(file=paste("figures/Comparison_Distance_",sp,"_",enh,".pdf",sep=""), width=6, height=6)
    par(mar=c(4.5, 4.5, 2.1, 1.1))

    dn=density(n$Distance)
    dc=density(c$min_dist_baitedTSS,bw=dn$bw)

    xlim=range(c(dn$x, dc$x))
    ylim=range(c(dn$y, dc$y))

    plot(dn$x, dn$y, xlim=xlim,ylim=ylim, xlab="distance", ylab="density", main="", type="l", col="black")
    lines(dc$x, dc$y, col="blue")

    legend("topright", inset=0.01, col=c("black", "blue"), legend=c("neighbors", "contacts"),lty=1)

    dev.off()

    ## proportion enhancers in common

    nb.common=unlist(lapply(common.genes, function(x) length(intersect(n$EnhancerID[which(n$GeneID==x)], c$enhancer[which(c$gene==x)]))))
    names(nb.common)=common.genes

    prop.common.contacts=100*nb.common[names(c.nbenh)]/c.nbenh
    prop.common.neighbors=100*nb.common[names(n.nbenh)]/n.nbenh

    pdf(file=paste("figures/CommonEnhancers_",sp,"_",enh,".pdf",sep=""), width=10, height=6)
    par(mfrow=c(1,2))
    hist(prop.common.contacts, xlab="% contacted enhancers found in neighbor regions", main="")
    hist(prop.common.neighbors, xlab="% neighboring enhancers found in contact", main="")
    dev.off()
    
  }
}

###############################################################################

