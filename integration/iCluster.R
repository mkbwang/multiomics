
rm(list=ls())


library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)

data(gbm)
mut.rate <- apply(gbm.mut, 2, mean)
gbm.mut2 <- gbm.mut[, which(mut.rate > 0.02)]

dim(gbm.exp)
data(variation.hg18.v10.nov.2010)


gbm.cn <- CNregions(seg=gbm.seg,epsilon=0,adaptive=FALSE,rmCNV=TRUE,
                 cnv=variation.hg18.v10.nov.2010[,3:5],
                 frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5)


gbm.cn=gbm.cn[order(rownames(gbm.cn)),]
all(rownames(gbm.cn)==rownames(gbm.exp))


fit.single <- iClusterPlus(dt1=gbm.mut2,dt2=gbm.cn,dt3=gbm.exp,
                        type=c("binomial","gaussian","gaussian"),
                        lambda=c(0.04,0.61,0.90),K=2,maxiter=10)



