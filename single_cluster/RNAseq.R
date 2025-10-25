
library(ggplot2)
library(patchwork)
rm(list=ls())


rna_normalized <- read.csv("data/RNA_normalized.csv", row.names=1) |> as.matrix()
rna_counts <- read.csv("data/RNA_counts.csv", row.names=1) |> as.matrix()

# visualize prevalences
prevalences <- data.frame(Prev=rowMeans(rna_counts > 0) * 100)
prev_plot <- ggplot(prevalences, aes(x=Prev)) +
    geom_histogram(breaks = seq(0, 100, by = 5), color = "black", fill = "white") +
    scale_x_continuous(limits = c(0, 100), breaks=seq(0, 100, 10))+
    xlab("Prevalence of Gene (%)") + ylab("Count") +
    ggtitle("RNASeq of 67 samples and 1362 genes") +
    theme_bw()


# focus on the genes whose prevalence is larger than the median
median_prevalence <- median(prevalences$Prev)
rna_normalized_subset <- rna_normalized[prevalences$Prev > median_prevalence, ] |> t()
rna_counts_subset <- rna_counts[prevalences$Prev > median_prevalence, ] |> t()


## dimension reduction
library(factoextra)
library(cluster)


rna_normalized_subset_scaled <- scale(rna_normalized_subset)

pca_result <- prcomp(rna_normalized_subset_scaled)
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

pca_variance_df <- data.frame(PC=seq(1, nrow(rna_normalized_subset_scaled)),
                              Variance=variance_explained*100)

pca_var_plot <- ggplot(pca_variance_df[1:20, ], aes(x=PC, y=Variance)) +
    geom_point(size=1.2) + geom_line() +
    scale_x_continuous(breaks=seq(1, 20), limits = c(1, 20))+
    scale_y_continuous(breaks=seq(0, 18, 3), limits=c(0, 16))+
    xlab("PC Number") + ylab("Variance Explained (%)")

pca_var_plot

## pick the top 3 PCs and do hierarchical clustering


top_scores <- pca_result$x[, 1:2]
sample_dist_mat <- dist(top_scores)
hclust_result <- hclust(sample_dist_mat, method="ward.D2")
colvecs <- c("#7DA1C4", "#9BB890", "#C59B75", "#9179AF",
             "#53ADA8", "#C5778D")
dgram <- fviz_dend(hclust_result, cex = 0.5, k=6,
          k_colors=c("#7DA1C4", "#9BB890", "#C59B75", "#9179AF",
                     "#53ADA8", "#C5778D"),
          main="Hierarchical Clustering Based on RNASeq",
          xlab="Sample", ylab="Euclidean Distance based on Two PCs")
clusters <- cutree(hclust_result, k = 6)


## visualze the PCA scores colored by cluster type
top_scores_df <- as.data.frame(top_scores)
top_scores_df$Cluster <- clusters
colvecs_reordered <- colvecs[top_scores_df[c("UPMC.20", "UPMC.11", "MS17.8290", "L45", "SU16.51836", "L64"), "Cluster"]]
top_scores_df$Cluster <- as.factor(top_scores_df$Cluster)

pc_plot <- ggplot(top_scores_df, aes(x=PC1, y=PC2, color=Cluster)) +
    geom_point() + scale_color_manual(values=colvecs_reordered)+
    xlab("PC1") + ylab("PC2")

write.csv(top_scores_df, "single_cluster/RNAseq_cluster.csv")
