rm(list=ls())
library(ggplot2)
library(patchwork)
library(multiviewtest)
library(vegan)
library(ComplexHeatmap)
library(mclust)
library(blockmodels)


# load RNASeq
rna_normalized <- read.csv("data/RNA_normalized.csv", row.names=1) |> as.matrix()
rna_counts <- read.csv("data/RNA_counts.csv", row.names=1) |> as.matrix()
rna_samples <- colnames(rna_normalized)


# load WES
wes <- read.csv("data/WES_data.csv", row.names=1)
wes_samples <- colnames(wes)


# subset samples
subset_samples <- intersect(rna_samples, wes_samples)
rna_counts_subset <- rna_counts[, subset_samples] |> as.matrix() |> t()
prevalences <- colMeans(rna_counts_subset > 0)
rna_normalized_subset <- rna_normalized[, subset_samples] |> as.matrix() |> t()
rna_normalized_subset <- rna_normalized_subset[, prevalences > 0.5]
rna_normalized_subset_scaled <- scale(rna_normalized_subset)
pca_result <- prcomp(rna_normalized_subset_scaled)
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
rna_pc_scores <- pca_result$x[, c(1,2)]


wes_subset <- wes[, subset_samples] |> as.matrix() |> t()
mutation_rate <- colMeans(wes_subset)
wes_subset <- wes_subset[, mutation_rate > 0.1]
wes_dist_mat <- dist(wes_subset, method="binary") |> as.matrix()

# visualize the jaccard distances
all_dists <- wes_dist_mat[upper.tri(wes_dist_mat)]
breaks <- seq(0, 1, 0.1)
binned <- cut(all_dists, breaks = breaks, include.lowest = TRUE)
df <- data.frame(bin = binned)
ggplot(df, aes(x = bin)) +
    geom_bar() +
    xlab("Jaccard Distance Range") +
    ylab("Count") +
    theme_bw()
wes_adjmat <- (wes_dist_mat < 0.8) * 1



# clustering
mclust_result_rna <- Mclust(data=rna_pc_scores, modelNames="VII")
optimal_rna_number <- length(unique(mclust_result_rna$classification))

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

combined_plots <- wrap_plots(rna_plot_list, ncol=2)

# stochastic block model
wes_sbm <- BM_bernoulli(membership_type="SBM", adj=wes_adjmat,
                        plotting="")
wes_sbm$estimate()
optimal_wes_number <- which.max(wes_sbm$ICL)
possible_wes_numbers <- seq(2, optimal_wes_number, 1)
wes_assignments <- matrix(0, nrow=nrow(wes_adjmat),
                          ncol=length(possible_wes_numbers))
for (j in possible_wes_numbers){

    cluster_result <- apply(wes_sbm$memberships[[j]]$Z, 1, function(vals){
        which.max(vals)
    })
    wes_assignments[, j-1] <- cluster_result

}
colnames(wes_assignments) <- sprintf("WES_%d", possible_wes_numbers)



# visualize binary adjacency matrices
for (j in possible_wes_numbers){

    cluster_assignments <- wes_assignments[, j-1]
    label_counts <- table(cluster_assignments)
    cum_counts <- cumsum(label_counts)
    reordered_adjmat <- wes_adjmat[order(cluster_assignments),
                                   order(cluster_assignments)]
    ht <- Heatmap(
        reordered_adjmat,
        name = "Adjacency",
        col = c("white", "black"),
        rect_gp = gpar(col = "grey80", lwd = 0.5),     # Light cell grid
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_heatmap_legend = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE
    )
    draw(ht)
    decorate_heatmap_body("Adjacency", {
        for (cut in cum_counts[-length(cum_counts)]) {
            grid.lines(x = c(0,1), y = rep(1 - cut/nrow(reordered_adjmat),2), gp = gpar(col="red", lwd=2))
        }
        for (cut in cum_counts[-length(cum_counts)]) {
            grid.lines(x = rep(cut/ncol(reordered_adjmat),2), y = c(0,1), gp = gpar(col="red", lwd=2))
        }
    })

}


cluster_df <- data.frame(Samples=rownames(rna_pc_scores),
                         RNA_score1=rna_pc_scores[, 1],
                         RNA_score2=rna_pc_scores[, 2])

cluster_df <- cbind(cluster_df, rna_assignments, wes_assignments)
write.csv(cluster_df, "multi_cluster/RNA_WES/RNA_WES_cluster.csv",
          row.names=FALSE)

