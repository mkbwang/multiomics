


rm(list=ls())


library(mclust)

rnaseq_cluster <- read.csv("single_cluster/RNAseq_cluster.csv")
TME_cluster <- read.csv("single_cluster/TME_cluster.csv")


wes_cluster <- read.csv("single_cluster/WES_cluster.csv")

overlap1 <- rnaseq_cluster %>% inner_join(TME_cluster, by="X")
ari_score <- adjustedRandIndex(overlap1$Cluster.x, overlap1$Cluster.y)

overlap2 <- rnaseq_cluster %>% inner_join(wes_cluster, by="X")
ari_score <- adjustedRandIndex(overlap2$Cluster.x, overlap2$Cluster.y)

overlap3 <- TME_cluster %>% inner_join(wes_cluster, by="X")
ari_score <- adjustedRandIndex(overlap3$Cluster.x, overlap3$Cluster.y)


IHC_cluster <- read.csv("single_cluster/IHC_cluster.csv")

