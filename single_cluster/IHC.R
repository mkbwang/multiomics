
rm(list=ls())
library(logisticPCA)
library(diceR)

IHC_data <- read.csv("data/IHC_percent.csv", row.names=1) |> na.omit() |> as.matrix()
IHC_binary <- read.csv("data/IHC_binary.csv", row.names=1) |> na.omit() |> as.matrix()

binary_trends <- unique(IHC_binary)
trend_groups <- list()


for (j in 1:nrow(binary_trends)){

    current_trend <- binary_trends[j, ]
    mask <- apply(IHC_binary, 1, identical, current_trend)
    trend_groups[[j]] <- cbind(IHC_binary[mask, ], j)

}

cluster_result <- do.call(rbind, trend_groups)

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
    breaks = c(0, 1),
    colors = c("white", "#333333")
)

ht <- Heatmap(cluster_result[, seq(1,4)], name="IHC", cluster_rows=FALSE, cluster_columns=FALSE,
              show_column_names = TRUE, show_row_names =  FALSE,
              show_heatmap_legend = FALSE,
              top_annotation = NULL, # Add column annotation
              left_annotation = NULL, # Add row annotation
              border_gp = gpar(col="black", lwd=0.5),
              width = unit(0.95, "npc"),
              height = unit(0.85, "npc"),
              na_col = "#555555",
              col = col_fun)

write.csv(as.data.frame(cluster_result),
          "single_cluster/IHC_cluster.csv", quote=FALSE)

