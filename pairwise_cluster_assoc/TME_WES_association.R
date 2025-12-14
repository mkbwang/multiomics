
rm(list=ls())
library(ggplot2)
library(patchwork)
library(multiviewtest)
library(vegan)
library(blockmodels)



# load TME
TME <- read.csv("data/TME.csv", row.names=1)
TME_neigh <- TME[, seq(1, 10)]
TME_neigh <- TME_neigh / 100
bray_dist <- vegdist(TME_neigh, method="bray") |> as.matrix()
TME_samples <- rownames(TME_neigh)

# load WES
wes <- read.csv("data/WES_data.csv", row.names=1)
wes_samples <- colnames(wes)


# subset samples
subset_samples <- intersect(wes_samples, TME_samples)
bray_dist_subset <- bray_dist[subset_samples, subset_samples]
mds_result <- cmdscale(d=bray_dist_subset, k=2,
                       eig=TRUE)
mds_coords <- mds_result$points


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
    xlab("Bray Curtis Distance Range") +
    ylab("Count") +
    theme_bw()

wes_adjmat <- (wes_dist_mat < 0.8) * 1


# clustering
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


cluster_df <- data.frame(Samples=rownames(mds_coords),
                         TME_score1=mds_coords[, 1],
                         TME_score2=mds_coords[, 2])

cluster_df <- cbind(cluster_df, TME_assignments, wes_assignments)

write.csv(cluster_df, "multi_cluster/TME_WES/TME_WES_cluster.csv",
          row.names=FALSE, quote=FALSE)


