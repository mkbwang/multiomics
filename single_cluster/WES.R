

library(ggplot2)
library(dplyr)
library(tidyr)
rm(list=ls())


wes <- read.csv("data/WES_data.csv", row.names=1)

gene_names <- rownames(wes)

# mutation rates
mutation_rate_df <- data.frame(Gene=rownames(wes),
                          Prevalence=rowSums(wes > 0))
rate_histogram <- ggplot(mutation_rate_df, aes(x=Prevalence)) +
    geom_histogram(breaks = seq(0, 40, by = 5), color = "black", fill = "white") +
    scale_x_continuous(limits = c(0, 36), breaks=seq(0, 35, 5))+
    scale_y_continuous(breaks=seq(0, 100, 5))+
    xlab("Number of Mutations across all 90 Samples") + ylab("Gene Count")+
    theme_bw()
## select features whose mutation rate is at least 10
wes_subset <- wes[mutation_rate_df$Prevalence >= 10, ] |> as.matrix()



# calculate feature correlation by fitting chi square tests
pvals_mat <- matrix(0, nrow=nrow(wes_subset), ncol=nrow(wes_subset))
rownames(pvals_mat) <- colnames(pvals_mat) <- rownames(wes_subset)
for (j in 1:(nrow(wes_subset)-1)){
    for (k in (j+1):nrow(wes_subset)){
        contingency_table <- table(wes_subset[j, ], wes_subset[k, ])
        test_result <- chisq.test(contingency_table)
        pvals_mat[j, k] <- pvals_mat[k, j] <- test_result$p.value
    }
}

## p value adjustment
padj_mat <- matrix(0, nrow=nrow(wes_subset), ncol=nrow(wes_subset))
rownames(padj_mat) <- colnames(padj_mat) <- rownames(wes_subset)

all_pvals <- pvals_mat[upper.tri(pvals_mat)]
all_padj <- p.adjust(all_pvals, method="BH")

padj_mat[upper.tri(padj_mat)] <- all_padj
padj_mat <- padj_mat + t(padj_mat)
diag(padj_mat) <- 1


significance_count <- colSums(padj_mat < 0.05)
## retain genes which have significant correlation with at least one other gene
padj_mat <- padj_mat[significance_count > 0, significance_count > 0]
library(igraph)
library(ggraph)
padj_graph <- graph_from_adjacency_matrix(
    1*(padj_mat < 0.05),
    mode = "undirected"
)
padj_graph_plot <- ggraph(padj_graph, layout = "kk") +  # 'kk' is a popular layout algorithm
    geom_edge_link(color = "gray") +
    geom_node_point(color = "salmon", size = 5) +
    geom_node_text(aes(label = name), vjust = 1.8) +
    theme_void() # Clean, empty theme


degree_df <- data.frame(Gene=names(significance_count),
                        Degree=significance_count)
mutation_info <- degree_df %>% left_join(mutation_rate_df, by="Gene")
write.csv(mutation_info, "single_cluster/Mutation_info.csv",
          row.names=FALSE, quote=FALSE)

selected_genes <- names(significance_count[significance_count > 0])
wes_subset <- wes_subset[selected_genes, ]



## logisticSVD, then hierarchical clustering

wes_subset <- t(wes_subset)
mutation_sum <- rowSums(wes_subset)

samples_nomutation <- rownames(wes_subset)[mutation_sum == 0]
wes_subset <- wes_subset[mutation_sum > 0, ]


library(logisticPCA)

lsvd_cv <- cv.lsvd(wes_subset, ks=seq(1, 5),
                   max_iters=5000)

cv_neglik <- data.frame(Rank=seq(1, 5),
                        NegLik=lsvd_cv[, 1])
cv_plot <- ggplot(cv_neglik, aes(x=Rank, y=NegLik)) + geom_point(size=1.2)+
    geom_line() + xlab("Rank of LogisticSVD") + ylab("Cross Validation Negative Log Likelihood")

lsvd_result <- logisticSVD(x=wes_subset, k=3)
scores <- lsvd_result$A
rownames(scores) <- rownames(wes_subset)

sample_dist_mat <- dist(scores)
hclust_result <- hclust(sample_dist_mat, method="ward.D2")
colvecs <- c("#7DA1C4", "#9BB890", "#C59B75", "#9179AF",
             "#53ADA8", "#C5778D")
dgram <- fviz_dend(hclust_result, cex = 0.5,k=6,
          k_colors=colvecs,
          main="Hierarchical Clustering Based on WES",
          xlab="Sample", ylab="Euclidean Distance based on LogisticSVD Scores")
clusters <- cutree(hclust_result, k = 6)

scores_df <- as.data.frame(scores)
scores_df$Cluster <- clusters
scores_df <- scores_df %>% arrange(Cluster)
scores_df$Color <- ""
scores_df["X17.18488", "Color"] <- "#7DA1C4"
scores_df["UPMC.13", "Color"] <- "#9BB890"
scores_df["UPMC.21", "Color"] <- "#C59B75"
scores_df["X18.797", "Color"] <- "#9179AF"
scores_df["X15.17492", "Color"] <- "#53ADA8"
scores_df["S12.7154", "Color"] <- "#C5778D"

color_cluster_match <- scores_df[c("X17.18488", "UPMC.13", "UPMC.21",
                                   "X18.797", "X15.17492", "S12.7154"), c("Cluster", "Color")] %>%
    arrange(Cluster)
reordered_color <- color_cluster_match$Color
scores_df$Color <- reordered_color[scores_df$Cluster]

null_df <- matrix(0, nrow=length(samples_nomutation), ncol=3) |> as.data.frame()
rownames(null_df) <- samples_nomutation
null_df$Cluster <- 7
null_df$Color <- "#222222"
final_result <- rbind(scores_df, null_df)
write.csv(final_result, "single_cluster/WES_cluster.csv")


#TODO: visualize heatmap
library(ComplexHeatmap)
library(circlize)
wes_subset_viz <- wes[colnames(wes_subset), ] |> t()
wes_subset_viz <- wes_subset_viz[rownames(final_result), ]



annot_df <- data.frame(
    Cluster = sprintf("Cluster%d", final_result$Cluster)
)
rownames(annot_df) <- rownames(wes_subset_viz)

pheno_colors <- unique(final_result$Color) |> unique()
names(pheno_colors) <- sprintf("Cluster%d", 1:7)
color_list <- list(
    Cluster = pheno_colors
)
ra <- rowAnnotation(
    df = annot_df,
    col = color_list,
    show_annotation_name = FALSE
)


col_fun <- colorRamp2(
    breaks = c(0, 1),
    colors = c("white", "#333333")
)

ht <- Heatmap(as.matrix(wes_subset_viz), name="Mutation", cluster_rows=FALSE, cluster_columns=FALSE,
              show_column_names = TRUE, show_row_names =  FALSE,
              show_heatmap_legend = FALSE,
              top_annotation = NULL, # Add column annotation
              left_annotation = ra, # Add row annotation
              border_gp = gpar(col="black", lwd=0.5),
              width = unit(0.95, "npc"),
              height = unit(0.85, "npc"),
              na_col = "#555555",
              col = col_fun)


