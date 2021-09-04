################################################################################

path="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
path1=paste(path, "backup/supplementary_datasets_backup_01092021/SupplementaryDataset5/",sep="")
path2=paste(path, "SupplementaryDataset5/",sep="")
pathRepeats="/beegfs/data/necsulea/RegulatoryLandscapes/results/sequence_composition/"

################################################################################

for(sp in c("human", "mouse")){
  rep=read.table(paste(pathRepeats, sp, "/restriction_fragments/overlap_exons_repeats.txt", sep=""), h=T, stringsAsFactors=F)
  rownames(rep)=rep$ID
  
  old.obs=read.table(paste(path1, sp, "/statistics_contacted_sequence_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  old.sim=read.table(paste(path1, sp, "/statistics_contacted_sequence_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  
  rownames(old.obs)=paste(old.obs$chr, old.obs$start, old.obs$end, sep=":")
  rownames(old.sim)=paste(old.sim$chr, old.sim$start, old.sim$end, sep=":")
  
  columns=intersect(c("chr", "start", "end", "length", "GC_content", "nb_baits_500kb", "nb_genes_500kb", "BLAT_match", "FANTOM5_bp", "ENCODE_bp", "RoadmapEpigenomics_bp", "FOCS_GRO_seq_bp"), colnames(old.obs))
  
  obs=old.obs[,columns]
  sim=old.sim[,columns]
  
  obs$repeat_bp=rep[rownames(old.obs), "RepeatLength"]
  obs$exon_bp=rep[rownames(old.obs), "ExonicLength"]
  obs$length=rep[rownames(old.obs), "TotalLength"]
  
  sim$repeat_bp=rep[rownames(old.sim), "RepeatLength"]
  sim$exon_bp=rep[rownames(old.sim), "ExonicLength"]
  sim$length=rep[rownames(old.sim), "TotalLength"]
  
  write.table(obs, file=paste(path2,  sp, "/statistics_contacted_sequence_original.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  write.table(sim, file=paste(path2,  sp, "/statistics_contacted_sequence_simulated.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)  
}

################################################################################
