
rm(list=ls())
library(ggplot2)
library(patchwork)
library(multiviewtest)
library(vegan)
library(blockmodels)


# load WES
wes <- read.csv("data/WES_data.csv", row.names=1)
wes_samples <- colnames(wes)


# load IHC
IHC_binary <- read.csv("data/IHC_binary.csv", row.names=1) |> na.omit() |> as.matrix()
IHC_samples <- rownames(IHC_binary)


# subset samples
subset_samples <- intersect(wes_samples, IHC_samples)
wes_subset <- wes[, subset_samples] |> as.matrix() |> t()
mutation_rate <- colMeans(wes_subset)
wes_subset <- wes_subset[, mutation_rate > 0.1]
wes_dist_mat <- dist(wes_subset, method="binary") |> as.matrix()
wes_adjmat <- (wes_dist_mat < 0.8) * 1
all_dists <- wes_dist_mat[upper.tri(wes_dist_mat)]
breaks <- seq(0, 1, 0.1)
binned <- cut(all_dists, breaks = breaks, include.lowest = TRUE)
df <- data.frame(bin = binned)
ggplot(df, aes(x = bin)) +
    geom_bar() +
    xlab("Bray Curtis Distance Range") +
    ylab("Count") +
    theme_bw()


IHC_subset <- IHC_binary[subset_samples, ]



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



cluster_df <- data.frame(Samples=rownames(IHC_subset))
cluster_df <- cbind(cluster_df, wes_assignments, IHC_subset)


write.csv(cluster_df, "multi_cluster/WES_IHC/WES_IHC_cluster.csv",
          row.names=FALSE, quote=FALSE)





