

rm(list=ls())
library(ggplot2)
library(patchwork)
library(multiviewtest)
library(vegan)
library(ComplexHeatmap)
library(mclust)
library(sbm)


# load RNASeq

rna_normalized <- read.csv("data/RNA_normalized.csv", row.names=1) |> as.matrix()
rna_counts <- read.csv("data/RNA_counts.csv", row.names=1) |> as.matrix()
rna_samples <- colnames(rna_normalized)

# load IHC
IHC_binary <- read.csv("data/IHC_binary.csv", row.names=1) |> na.omit() |> as.matrix()
IHC_samples <- rownames(IHC_binary)



# subset samples
subset_samples <- intersect(rna_samples, IHC_samples)
rna_counts_subset <- rna_counts[, subset_samples] |> as.matrix() |> t()
prevalences <- colMeans(rna_counts_subset > 0)
rna_normalized_subset <- rna_normalized[, subset_samples] |> as.matrix() |> t()
rna_normalized_subset <- rna_normalized_subset[, prevalences > 0.5]
rna_normalized_subset_scaled <- scale(rna_normalized_subset)
pca_result <- prcomp(rna_normalized_subset_scaled)
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
rna_pc_scores <- pca_result$x[, c(1,2)]

IHC_subset <- IHC_binary[subset_samples, ]



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




cluster_df <- data.frame(Samples=rownames(rna_pc_scores),
                         RNA_score1=rna_pc_scores[, 1],
                         RNA_score2=rna_pc_scores[, 2])

cluster_df <- cbind(cluster_df, rna_assignments, IHC_subset)

write.csv(cluster_df, "multi_cluster/RNA_IHC/RNA_IHC_cluster.csv",
          row.names=FALSE, quote=FALSE)

