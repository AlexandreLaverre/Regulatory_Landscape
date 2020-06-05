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
shhchr="7"

###########################################################################

baitcoords=read.table(paste(pathInteractions, tolower(sp), "/bait_coords_",assembly,".txt",sep=""), h=T, stringsAsFactors=F)
shhbait=baitcoords[grep(shhid, baitcoords$gene_ID), ]
  
###########################################################################

interactions=read.table(paste(pathInteractions, tolower(sp), "/all_interactions.txt",sep=""), h=T, stringsAsFactors=F)

## select only cis interactions, unbaited, within the distance range

interactions=interactions[which(interactions$chr_bait==interactions$chr),]
interactions=interactions[which(interactions$type=="unbaited"),]
interactions=interactions[which(interactions$distance>=minDistance & interactions$distance<=maxDistance),]

shhinteractions=interactions[which(interactions$chr_bait==shhbait$chr & interactions$start_bait==shhbait$start & interactions$end_bait==shhbait$end),]
  
###########################################################################

## select limits for the plot

margin=10000

shhxlim=range(c(shhinteractions$start, shhinteractions$end, shhbait$start, shhbait$end))+c(-margin, margin)

shhgenecoords=genecoords[which(genecoords$chr==shhchr & ((genecoords$start>=shhxlim[1] & genecoords$start<=shhxlim[2]) | (genecoords$end>=shhxlim[1] & genecoords$end<=shhxlim[2]) |  (genecoords$start<=shhxlim[1] & genecoords$end>=shhxlim[2]))),]

## other baits in regions

allshhbaits=baitcoords[which(baitcoords$chr==shhchr & baitcoords$start>=shhxlim[1] & baitcoords$end<=shhxlim[2]),]

###########################################################################

save(list=c("allshhbaits", "shhid", "shhchr", "shhinteractions", "shhxlim", "shhgenecoords"), file="RData/data.Shh.figure.RData")

###########################################################################
