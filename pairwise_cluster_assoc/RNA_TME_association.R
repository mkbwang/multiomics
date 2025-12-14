

rm(list=ls())
library(ggplot2)
library(patchwork)
library(multiviewtest)
library(vegan)
library(mclust)

# load RNASeq

rna_normalized <- read.csv("data/RNA_normalized.csv", row.names=1) |> as.matrix()
rna_counts <- read.csv("data/RNA_counts.csv", row.names=1) |> as.matrix()
rna_samples <- colnames(rna_normalized)


# load TME
TME <- read.csv("data/TME.csv", row.names=1)
TME_neigh <- TME[, seq(1, 10)]
TME_neigh <- TME_neigh / 100
bray_dist <- vegdist(TME_neigh, method="bray") |> as.matrix()
TME_samples <- rownames(TME_neigh)


# subset samples
subset_samples <- intersect(rna_samples, TME_samples)
rna_counts_subset <- rna_counts[, subset_samples] |> as.matrix() |> t()
prevalences <- colMeans(rna_counts_subset > 0)
rna_normalized_subset <- rna_normalized[, subset_samples] |> as.matrix() |> t()
rna_normalized_subset <- rna_normalized_subset[, prevalences > 0.5]
rna_normalized_subset_scaled <- scale(rna_normalized_subset)
pca_result <- prcomp(rna_normalized_subset_scaled)
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
rna_pc_scores <- pca_result$x[, c(1,2)]


bray_dist_subset <- bray_dist[subset_samples, subset_samples]
## multidimensional scaling
mds_result <- cmdscale(d=bray_dist_subset, k=2,
                       eig=TRUE)
mds_coords <- mds_result$points



# clustering
mclust_result_rna <- Mclust(data=rna_pc_scores, modelNames="VII")
optimal_rna_number <- length(unique(mclust_result_rna$classification))

mclust_result_TME <- Mclust(data=mds_coords, modelNames="VII")
optimal_TME_number <- length(unique(mclust_result_TME$classification))


possible_rna_numbers <- seq(2, optimal_rna_number, 1)
rna_assignments <- matrix(0, nrow=nrow(rna_pc_scores),
                          ncol=length(possible_rna_numbers))
rna_plot_list <- list()
for (j in possible_rna_numbers){
    cluster_result <- Mclust(data=rna_pc_scores, modelNames="VII", G=j)
    rna_assignments[, j-1] <- cluster_result$classification
    rna_df <- data.frame(RNA_cluster=as.factor(cluster_result$classification),
                         RNA_score1=rna_pc_scores[, 1],
                         RNA_score2=rna_pc_scores[, 2])
    rna_plot_list[[j-1]] <- ggplot(rna_df, aes(x=RNA_score1, y=RNA_score2, color=RNA_cluster)) +
        geom_point(alpha=0.7, size=2) + scale_color_manual(values=c("#004949", "#490092", "#920000", "#924900",
                                                                    "#333333")) +
        theme_bw() + theme(legend.position = "none") +
        xlab("PC1") + ylab("PC2") + ggtitle(sprintf("%d Clusters", j))
}
colnames(rna_assignments) <- sprintf("RNA_%d", possible_rna_numbers)


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


cluster_df <- data.frame(Samples=rownames(rna_pc_scores),
                         RNA_score1=rna_pc_scores[, 1],
                         RNA_score2=rna_pc_scores[, 2],
                         TME_score1=mds_coords[, 1],
                         TME_score2=mds_coords[, 2])
cluster_df <- cbind(cluster_df, rna_assignments, TME_assignments)
write.csv(cluster_df, "multi_cluster/RNA_TME/RNA_TME_cluster.csv",
          row.names=FALSE)


combined_plots <- wrap_plots(TME_plot_list, ncol=2)


