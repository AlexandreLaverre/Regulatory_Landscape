#############################" Corrélation entre les variables ############################
library("corrplot")
conserv <- read.table("PIR_cons_all_overlap.txt", header=T)
simul <- read.table("PIR_cons_all_overlap_simul.txt", header=T)

conserv2 <- conserv[,-1]
rownames(conserv2) <- conserv[,1]
M <- cor(conserv2, method="pearson")
cor.mtest <- function(conserv2, method="pearson") {
  conserv2 <- as.matrix(conserv2)
  n <- ncol(conserv2)
  p.conserv2<- matrix(NA, n, n)
  diag(p.conserv2) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(conserv2[, i], conserv2[, j], method="pearson")
      p.conserv2[i, j] <- p.conserv2[j, i] <- tmp$p.value
    }
  }
  colnames(p.conserv2) <- rownames(p.conserv2) <- colnames(conserv2)
  p.conserv2
}
# Matrice de p-value de la corrélation
p.conserv2 <- cor.mtest(conserv2)
head(p.conserv2[, 1:5])

corrplot(M, method="circle", type="upper", tl.col="black", order = "hclust", tl.srt=45, p.mat = p.conserv2, sig.level = 0.01)
corrplot(M, method="circle", type="upper", tl.col="black", tl.srt=45, p.mat = p.conserv2, sig.level = 0.01)

#pairs(conserv2, cex=0.1)