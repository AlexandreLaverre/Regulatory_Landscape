################################################################################

source("../main_figures/parameters.R") ## paths are defined based on the user name

path1=paste(pathFinalData, "backup/supplementary_datasets_backup_01092021/SupplementaryDataset4/",sep="")
path2=paste(pathFinalData, "/tmp_writtable/SupplementaryDataset4/",sep="")

################################################################################

enh.datasets=list()
enh.datasets[["human"]]=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")
enh.datasets[["mouse"]]=c("ENCODE", "FANTOM5")

################################################################################

for(sp in c("human", "mouse")){
  for(enh in enh.datasets[[sp]]){
    print(sp, enh)
    rep=read.table(paste(pathRepeats, sp, "/", enh, "/overlap_exons_repeats.txt", sep=""), h=T, stringsAsFactors=F)
    rownames(rep)=rep$ID

    old.obs=read.table(paste(path1, sp, "/", enh, "/statistics_contacted_enhancers_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    old.sim=read.table(paste(path1, sp, "/", enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    rownames(old.obs)=paste(old.obs$chr, old.obs$start, old.obs$end, sep=":")
    rownames(old.sim)=paste(old.sim$chr, old.sim$start, old.sim$end, sep=":")

    obs=old.obs[,c("chr", "start", "end", "length", "median_dist", "GC_content", "nb_baits_500kb", "nb_genes_500kb", "BLAT_match")]
    sim=old.sim[,c("chr", "start", "end", "length", "median_dist", "GC_content", "nb_baits_500kb", "nb_genes_500kb", "BLAT_match")]

    obs$repeat_bp=rep[rownames(old.obs), "RepeatLength"]
    obs$exon_bp=rep[rownames(old.obs), "ExonicLength"]

    sim$repeat_bp=rep[rownames(old.sim), "RepeatLength"]
    sim$exon_bp=rep[rownames(old.sim), "ExonicLength"]

    write.table(obs, file=paste(path2,  sp, "/", enh, "/statistics_contacted_enhancers_original.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
    write.table(sim, file=paste(path2,  sp, "/", enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)  
    
  }
}

################################################################################
