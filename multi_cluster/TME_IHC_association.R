

rm(list=ls())
library(ggplot2)
library(patchwork)
library(multiviewtest)
library(vegan)



# load TME
TME <- read.csv("data/TME.csv", row.names=1)
TME_neigh <- TME[, seq(1, 10)]
TME_neigh <- TME_neigh / 100
bray_dist <- vegdist(TME_neigh, method="bray") |> as.matrix()
TME_samples <- rownames(TME_neigh)


# load IHC
IHC_binary <- read.csv("data/IHC_binary.csv", row.names=1) |> na.omit() |> as.matrix()
IHC_samples <- rownames(IHC_binary)



# subset samples
subset_samples <- intersect(TME_samples, IHC_samples)
bray_dist_subset <- bray_dist[subset_samples, subset_samples]
mds_result <- cmdscale(d=bray_dist_subset, k=2,
                       eig=TRUE)
mds_coords <- mds_result$points
IHC_subset <- IHC_binary[subset_samples, ]


# GMM clustering
mclust_result_TME <- Mclust(data=mds_coords, modelNames="VII")
optimal_TME_number <- length(unique(mclust_result_TME$classification))

possible_TME_numbers <- seq(2, optimal_TME_number, 1)
TME_assignments <- matrix(0, nrow=nrow(mds_coords),
                          ncol=length(possible_TME_numbers))
TME_plot_list <- list()
for (j in possible_TME_numbers){

    cluster_result <- Mclust(data=mds_coords, modelNames="VII", G=j)
    TME_assignments[, j-1] <- cluster_result$classification
    TME_df <- data.frame(TME_cluster=as.factor(cluster_result$classification),
                         TME_score1=mds_coords[, 1],
                         TME_score2=mds_coords[, 2])
    TME_plot_list[[j-1]] <- ggplot(TME_df, aes(x=TME_score1, y=TME_score2, color=TME_cluster)) +
        geom_point(alpha=0.7, size=2) + scale_color_manual(values=c("#004949", "#490092", "#920000", "#924900",
                                                                    "#333333")) +
        theme_bw() + theme(legend.position = "none") +
        xlab("PC1") + ylab("PC2") + ggtitle(sprintf("%d Clusters", j))

}
colnames(TME_assignments) <- sprintf("TME_%d", possible_TME_numbers)

combined_plots <- wrap_plots(TME_plot_list, ncol=2)





cluster_df <- data.frame(Samples=rownames(mds_coords),
                         TME_score1=mds_coords[, 1],
                         TME_score2=mds_coords[, 2])

cluster_df <- cbind(cluster_df, TME_assignments, IHC_subset)



write.csv(cluster_df, "multi_cluster/TME_IHC/TME_IHC_cluster.csv",
          row.names=FALSE, quote=FALSE)

