
library(ggplot2)
library(tidyr)
rm(list=ls())
TME <- read.csv("data/TME.csv", row.names=1)


TME_neigh <- TME[, seq(1, 10)]
TME_neigh_long <- pivot_longer(TME_neigh, cols=colnames(TME_neigh),
                               names_to="Neighborhood",
                               values_to="Proportion")

composition_plot <- ggplot(TME_neigh_long, aes(x=Neighborhood, y=Proportion)) +
    geom_boxplot(width=0.1, fill="white")+
    labs(x="Neighborhood", y="Proportion (%)")+
    scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 10))+
    theme_bw()

TME_neigh <- TME_neigh / 100


# calculate bray-curtis distance
library(vegan)
library(ape)

bray_dist <-  vegdist(TME_neigh, method="bray")

## multidimensional scaling
mds_result <- cmdscale(d=bray_dist, k=2,
                       eig=TRUE)
mds_coords <- mds_result$points



## hierarchical clustering
library(factoextra)
library(cluster)

hclust_result <- hclust(bray_dist, method="average")
colvecs <- c("#7DA1C4", "#9BB890", "#C59B75", "#9179AF",
             "#C5778D")
dgram <- fviz_dend(hclust_result, cex = 0.5,
          k=5,
          k_colors=colvecs,
          main="Hierarchical Clustering Based on TME Neighborhood Proportion",
          xlab="Sample", ylab="Bray Curtis Distance")

clusters <- cutree(hclust_result, k = 5)
mds_coords_df <- as.data.frame(mds_coords)
mds_coords_df$Cluster <- clusters

colvecs_reordered <- colvecs[mds_coords_df[c("S13.945", "L64", "UPMC.20", "UPMC.11", "UPMC.21"), "Cluster"]]
mds_coords_df$Cluster <- as.factor(mds_coords_df$Cluster)
pcoa_plot <- ggplot(mds_coords_df, aes(x=V1, y=V2, color=Cluster)) +
    geom_point() + scale_color_manual(values=colvecs_reordered)+
    xlab("PCoA1") + ylab("PCoA2")

write.csv(mds_coords_df, "single_cluster/TME_cluster.csv")
