
enh.closest <- read.table(paste(pathFinalData, "/ENCODE.original2closest.TSS.txt", sep=""), h=T)

enh.closest$dist_class <- cut(enh.closest$closest.enh, breaks=seq(from=0, to=1e6, by=20000), include.lowest = T)
mean.val=tapply(100*enh.closest$conserv_mouse, enh.closest$dist_class, function(x) mean(x, na.rm=T))
ci.low=tapply(100*enh.closest$conserv_mouse, enh.closest$dist_class, function(x) t.test(x)[["conf.int"]][1])
ci.high=tapply(100*enh.closest$conserv_mouse, enh.closest$dist_class, function(x) t.test(x)[["conf.int"]][2])

enh.closest.sim <- read.table(paste(pathFinalData, "/ENCODE.simulated2closest.TSS.txt", sep=""), h=T)

enh.closest.sim$dist_class <- cut(enh.closest.sim$closest.enh, breaks=seq(from=0, to=1e6, by=20000), include.lowest = T)
mean.val.sim=tapply(100*enh.closest.sim$conserv_mouse, enh.closest.sim$dist_class, function(x) mean(x, na.rm=T))
ci.low.sim=tapply(100*enh.closest.sim$conserv_mouse, enh.closest.sim$dist_class, function(x) t.test(x)[["conf.int"]][1])
ci.high.sim=tapply(100*enh.closest.sim$conserv_mouse, enh.closest.sim$dist_class, function(x) t.test(x)[["conf.int"]][2])


ylim=c(0, 100)
dy=diff(ylim)/20
ylim=ylim+c(-dy, dy)

nbclasses=length(levels(enh.closest$dist_class))
xpos=1:nbclasses
xlim=c(0, max(xpos)+1)
class_leg <- c("0", "150", "300", "450", "600")
xax=seq(from=1, to=max(xpos), by=3)

plot(1, type="n", xlab="", ylab="", axes=F, main="ENCODE human", xlim=xlim, ylim=ylim, xaxs="i")

lines(xpos, mean.val, col=dataset.colors["Original"])
segments(xpos, ci.low, xpos, ci.high, col=dataset.colors["Original"])

lines(xpos, mean.val.sim, col=dataset.colors["Simulated"])
segments(xpos, ci.low.sim, xpos, ci.high.sim, col=dataset.colors["Simulated"])

axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)
mtext("distance to closest TSS (kb)", side=1, line=2.2)

axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
mtext("% aligned sequence", side=2, line=3)

##########################
load(paste(pathFigures, "RData/data.sequence.conservation.pcungapped.", ref_sp, ".Rdata", sep=""))

enh_align_obs$density_class <- cut(frag_align_obs$density50kb, breaks=seq(from=minDistance, to=maxDistance, by=1000), include.lowest = T)

nbclasses=length(levels(frag_align_obs$dist_class))
xpos=1:nbclasses

xlim=c(-0.5, max(xpos)+1)

## axis position
class_leg <- c("0", "0.5", "1", "1.5", "2")
xax=seq(from=0, to=max(xpos)+1, by=10)

labels=c("e", "g")
names(labels)=c(close_sp, target_sp)
other_sp = "mouse"

mean.val.obs=tapply(100*align_enhancer_obs[, other_sp], align_enhancer_obs$dist_class, function(x) mean(x, na.rm=T))
ci.low.obs=tapply(100*align_enhancer_obs[, other_sp], align_enhancer_obs$dist_class, function(x) t.test(x)[["conf.int"]][1])
ci.high.obs=tapply(100*align_enhancer_obs[, other_sp], align_enhancer_obs$dist_class, function(x) t.test(x)[["conf.int"]][2])

mean.val.simul=tapply(100*align_enhancer_simul[, other_sp], align_enhancer_simul$dist_class, function(x) mean(x, na.rm=T))
ci.low.simul=tapply(100*align_enhancer_simul[, other_sp], align_enhancer_simul$dist_class, function(x) t.test(x)[["conf.int"]][1])
ci.high.simul=tapply(100*align_enhancer_simul[, other_sp], align_enhancer_simul$dist_class, function(x) t.test(x)[["conf.int"]][2])

ylim=range(c(ci.low.obs, ci.high.obs, ci.low.simul, ci.high.simul))

dy=diff(ylim)/20
ylim=ylim+c(-dy, dy)

plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")

lines(xpos, mean.val.obs, col=dataset.colors["Original"])
segments(xpos, ci.low.obs, xpos, ci.high.obs, col=dataset.colors["Original"])

lines(xpos, mean.val.simul, col=dataset.colors["Simulated"])
segments(xpos, ci.low.simul, xpos, ci.high.simul, col=dataset.colors["Simulated"])

axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)
mtext("distance to promoters (Mb)", side=1, line=2.2, cex=0.8)

axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
mtext("% aligned sequence", side=2, line=3, cex=0.8)

mtext(paste(ref_sp, " vs. ", other_sp, ", enhancers", sep=""), side=3, cex=0.8)

mtext(labels[other_sp], side=3, line=1, at=-7.75, font=2, cex=1.2)


