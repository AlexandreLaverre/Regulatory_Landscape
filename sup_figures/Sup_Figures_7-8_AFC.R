library(ade4)
library(factoextra)
library(ggpubr)

ref_sp = "mouse"
path <- "/home/laverre/Data/Regulatory_landscape/result/"
original_data <- read.table(paste(path, "Supplementary_dataset1_original_interactions/", ref_sp, "/all_interactions.txt", sep=""), header = T)


#################### AFC ####################
##### Filters 
test_data <- original_data[which(original_data$type == "unbaited" & original_data$distance >= 25000 & original_data$distance <= 10000000),]

data <- test_data[,-c(1:8)]
data[!is.na(data)] <- 1
data[is.na(data)] <- 0

data_t <- data.frame(t(data)) # row = samples

AFC <- dudi.coa(data_t, scannf=F, nf=3)
AFC
summary(AFC)
AFC$eig
screeplot(AFC)

# Dendrogram
d=dist.dudi(AFC) 
cah=hclust(d,"ward")

dendo <- fviz_dend(cah, cex=0.5)
sample <- fviz_ca_row(AFC, repel = TRUE)

# Plot
pdf(paste(path, "Figures/AFC_", ref_sp, ".pdf", sep=""), width=8, height=4)
ggarrange(dendo, sample, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()


#plot(cah, main="Cluster Dendrogram (filtered)",hang=-1)

# Projection
#scatter(AFC, clab.row=0, clab.col=0.5, posieig="none")
#s.label(AFC$li, boxes=F, sub="Samples", ylim=c(-0.5,1.8), xlim=c(-2,1))

