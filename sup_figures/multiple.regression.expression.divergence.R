###############################################################################
library(car)

source("../main_figures/parameters.R")

ref_sp="mouse"

enh="ENCODE"

load(paste(pathFigures, "RData/data.", ref_sp, ".CM2019.AllOrgans.expdiv.RData", sep=""))
load(paste(pathFigures, "RData/data.regland.conservation.RData", sep=""))

regcons=regland.conservation[[ref_sp]][[enh]]

##############################################################################

common = intersect(rownames(expdiv), rownames(regcons))

gene.evolution=data.frame("ExpressionDivergence"=expdiv[common, "CorrelationSpearman"], stringsAsFactors=F, row.names=common)

gene.evolution$ExpressionLevel = expdiv[common, "MeanRPKM"]
gene.evolution$Specificity = expdiv[common, "TauHuman"]

gene.evolution$NbContacts = regcons[common, "class.nb.contacts.all"]
gene.evolution$MeanAlignment = regcons[common, "class.aln.score.all"]
#gene.evolution$MeanAlignment = regcons[common, "mean.phyloP.score.all"]
gene.evolution$Synteny = regcons[common, "class.synteny.cons.all"]
gene.evolution$ContactConservation = regcons[common, "class.contact.cons.all"]

##############################################################################

#### Model
mod=lm(ExpressionDivergence ~ ExpressionLevel + Specificity + NbContacts +
         MeanAlignment + Synteny + ContactConservation, data=gene.evolution)

summary(mod)
Anova(mod)
AIC(mod)
par(mfrow=c(2,2))
plot(mod)

#### Model with interactions
mod=lm(ExpressionDivergence ~ ExpressionLevel + Specificity + NbContacts +
         MeanAlignment * Synteny * ContactConservation, data=gene.evolution)

summary(mod)
AIC(mod)
par(mfrow=c(2,2))
plot(mod)
