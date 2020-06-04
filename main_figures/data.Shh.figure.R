###########################################################################

pathAnnotations="../../data/ensembl_annotations/"
pathInteractions="../../../RegulatoryLandscapesManuscript/SupplementaryDataset1/"

###########################################################################

sp="Human"

ensrelease=94
assembly="hg38"
minDistance=25e3
maxDistance=2.5e6

###########################################################################

genecoords=read.table(paste(pathAnnotations, sp, "/GeneInfo_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
colnames(genecoords)=c("id", "biotype", "description", "chr", "start", "end", "strand")

###########################################################################

shhid="ENSG00000164690" ## Shh
sshchr="7"

###########################################################################

baitcoords=read.table(paste(pathInteractions, tolower(sp), "/bait_coords_",assembly,".txt",sep=""), h=T, stringsAsFactors=F)
sshbait=baitcoords[grep(shhid, baitcoords$gene_ID), ]
  
###########################################################################

interactions=read.table(paste(pathInteractions, tolower(sp), "/all_interactions.txt",sep=""), h=T, stringsAsFactors=F)

## select only cis interactions, unbaited, within the distance range

interactions=interactions[which(interactions$chr_bait==interactions$chr),]
interactions=interactions[which(interactions$type=="unbaited"),]
interactions=interactions[which(interactions$distance>=minDistance & interactions$distance<=maxDistance),]

sshinteractions=interactions[which(interactions$chr_bait==sshbait$chr & interactions$start_bait==sshbait$start & interactions$end_bait==sshbait$end),]
  
###########################################################################

## select limits for the plot

sshxlim=range(c(sshinteractions$start, sshinteractions$end, sshbait$start, sshbait$end))

sshgenecoords=genecoords[which(genecoords$chr==sshchr & ((genecoords$start>=sshxlim[1] & genecoords$start<=sshxlim[2]) | (genecoords$end>=sshxlim[1] & genecoords$end<=sshxlim[2]) |  (genecoords$start<=sshxlim[1] & genecoords$end>=sshxlim[2]))),]

###########################################################################

save(list=c("baitcoords", "sshid", "sshchr", "sshinteractions", "sshxlim", "sshgenecoords"), file="RData/data.Shh.figure.RData")

###########################################################################
