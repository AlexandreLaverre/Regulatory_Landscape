options(stringsAsFactors = FALSE)

path <- "/home/laverre/Documents/Regulatory_Landscape/data/"
ref_sp = "mouse"

overlap <- read.table(paste(path, ref_sp, "/potential_enhancers/enhancers_coordinates/common_enhancers_overlap", sep=""), header=T, row.names=1)

overlap[!is.na(overlap)] <- TRUE
overlap[is.na(overlap)] <- FALSE

colnames(overlap) <- c("ENCODE", "FANTOM5")
overlap$ENCODE <- as.logical(overlap$ENCODE)
overlap$FANTOM5 <- as.logical(overlap$FANTOM5)

if (ref_sp == "human"){
  colnames(overlap) <- c("ENCODE", "FANTOM5", "RoadMap", "GRO-seq")
  overlap$RoadMapEpigenomics <- as.logical(overlap$RoadMap)
  overlap$GRO.seq <- as.logical(overlap$GRO-seq)
}


library(eulerr)
pdf(paste("/home/laverre/Data/Regulatory_landscape/result/Figures/Eulerr_diagram_enhancers_", ref_sp, ".pdf", sep=""))
plot(euler(overlap), quantities = TRUE, fill = c("skyblue", "tomato", "seagreen3", "tan1"))
dev.off()
