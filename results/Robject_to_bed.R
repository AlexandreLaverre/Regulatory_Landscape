setwd("/home/laverre/Téléchargements")
df <- data.frame(seqnames=seqnames(dreg.bs),
                 starts=start(dreg.bs)-1,
                 ends=end(dreg.bs),
                 names=c(rep(".", length(dreg.bs))),
                 scores=c(rep(".", length(dreg.bs))),
                 strands=strand(dreg.bs))

write.table(df, file="FANTOM5_enhancer_genomic_positions_hg19.bed", quote=F, sep="\t", row.names=F, col.names=F)
