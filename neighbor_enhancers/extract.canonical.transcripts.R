##########################################################################

pathEnsembl="../../data/ensembl_annotations/"
pathResults="../../results/neighbor_enhancers/"

release=94

##########################################################################

library(seqinr)

##########################################################################

for(sp in c("Human", "Mouse")){
  ## read transcript CDS sequences

  cds=read.fasta(paste(pathEnsembl, sp, "/AllCDS_Ensembl", release, ".fa",sep=""), seqtype="DNA", forceDNAtolower=FALSE)

  lencds=unlist(lapply(cds, length))
  names(lencds)=unlist(lapply(names(lencds), function(x) unlist(strsplit(x, split="\\."))[1]))

  ## read transcript cDNA sequences

  cdna=read.fasta(paste(pathEnsembl, sp, "/AllTranscripts_Ensembl",release,"_noMT.fa",sep=""), seqtype="DNA", forceDNAtolower=FALSE)
  lencdna=unlist(lapply(cdna, length))
  names(lencdna)=unlist(lapply(names(lencdna), function(x) unlist(strsplit(x, split=":"))[2]))

  ## APPRIS annotations
  appris=read.table(paste(pathEnsembl, sp, "/APPRIS_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(appris)=appris$Transcript.stable.ID

  ## transcript info
  txinfo=read.table(paste(pathEnsembl, sp, "/TranscriptInfo_Ensembl", release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  colnames(txinfo)[1]="gene_id"
  colnames(txinfo)[2]="transcript_id"
  colnames(txinfo)[4]="chr"
  colnames(txinfo)[5]="start"
  colnames(txinfo)[6]="end"
  colnames(txinfo)[7]="strand"
  
  ## gene info
  geneinfo=read.table(paste(pathEnsembl, sp, "/GeneInfo_Ensembl", release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(geneinfo)=geneinfo$stable_id

  ## add gene biotype to transcript info
  txinfo$gene_biotype=geneinfo[txinfo$gene_id,"biotype"]
  
  ## select protein-coding transcripts for protein-coding genes
  ## all transcripts for lncRNA genes
  ## we don't keep the other genes
  
  txinfo=txinfo[which((txinfo$biotype=="protein_coding" & txinfo$gene_biotype=="protein_coding") | (txinfo$gene_biotype%in%c("lincRNA", "antisense"))),]

  ## add transcript length

  txinfo$length=rep(NA, dim(txinfo)[1])
  txinfo$length[which(txinfo$biotype=="protein_coding")]=lencds[txinfo$transcript_id[which(txinfo$biotype=="protein_coding")]]
  txinfo$length[which(txinfo$biotype!="protein_coding")]=lencdna[txinfo$transcript_id[which(txinfo$biotype!="protein_coding")]]

  ## remove NA transcript lengths - haplotypes

  txinfo=txinfo[which(!is.na(txinfo$length)),]

  ## order transcripts by length

  txinfo=txinfo[order(txinfo$length, decreasing=T),]

  ## put "principal transcripts" first

  txinfo$APPRIS=appris[txinfo$transcript_id, "APPRIS.annotation"]
  p1=which(txinfo$APPRIS=="principal1")
  p2=which(txinfo$APPRIS=="principal2")
  p3=which(txinfo$APPRIS=="principal3")
  p4=which(txinfo$APPRIS=="principal4")
  p5=which(txinfo$APPRIS=="principal5")
  a1=which(txinfo$APPRIS=="alternative1")
  a2=which(txinfo$APPRIS=="alternative2")
  
  other=setdiff(1:dim(txinfo)[1], c(p1, p2, p3, p4, p5, a1, a2))
  txinfo=txinfo[c(p1, p2, p3, p4, p5, a1, a2, other),]

  ## select the first transcript for each gene
  ## APPRIS principal for the genes that have them
  ## longest CDS or longest cDNA for the other ones

  first=which(!duplicated(txinfo$gene_id))
  txinfo=txinfo[first,]

  ## extract TSS

  txinfo$TSS=rep(NA, dim(txinfo)[1])
  txinfo$TSS[which(txinfo$strand==1)]=txinfo$seq_region_start[which(txinfo$strand==1)]
  txinfo$TSS[which(txinfo$strand==-1)]=txinfo$seq_region_end[which(txinfo$strand==-1)]

  write.table(txinfo, file=paste(pathResults, tolower(sp), "canonical_transcripts_Ensembl",release,".txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  
}

##########################################################################
