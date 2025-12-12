
library(ggplot2)
library(patchwork)
library(tidyr)
library(ComplexHeatmap)
rm(list=ls())

## RNA
rna_normalized <- read.csv("data/RNA_normalized.csv", row.names=1) |> as.matrix()
rna_counts <- read.csv("data/RNA_counts.csv", row.names=1) |> as.matrix()


rna_prevalences <- data.frame(Prev=rowMeans(rna_counts > 0) * 100)
rna_prev_plot <- ggplot(rna_prevalences, aes(x=Prev)) +
    geom_histogram(breaks = seq(0, 100, by = 5), color = "black", fill = "white") +
    scale_x_continuous(limits = c(0, 100), breaks=seq(0, 100, 10))+
    xlab("Prevalence of Gene (%)") + ylab("Count") +
    ggtitle("RNASeq of 67 samples and 1362 genes") +
    theme_bw()
write.csv(rna_prevalences, "EDA_plots/RNAseq_prevalence.csv")


## WES
wes <- read.csv("data/WES_data.csv", row.names=1)
gene_names <- rownames(wes)

mutation_rate_df <- data.frame(Gene=rownames(wes),
                               Prevalence=rowSums(wes > 0))
mutation_rate_histogram <- ggplot(mutation_rate_df, aes(x=Prevalence)) +
    geom_histogram(breaks = seq(0, 40, by = 5), color = "black", fill = "white") +
    scale_x_continuous(limits = c(0, 36), breaks=seq(0, 35, 5))+
    scale_y_continuous(breaks=seq(0, 100, 5))+
    xlab("Number of Mutations across all 90 Samples") + ylab("Gene Count")+
    theme_bw()
write.csv(mutation_rate_df, "EDA_plots/wes_mutation.csv")


## TME
TME <- read.csv("data/TME.csv", row.names=1)

TME_neigh <- TME[, seq(1, 10)]
TME_neigh_long <- pivot_longer(TME_neigh, cols=colnames(TME_neigh),
                               names_to="Neighborhood",
                               values_to="Proportion")


TME_neigh_avg <- data.frame(Relabd=colMeans(TME_neigh),
                            Neighborhood=colnames(TME_neigh))
TME_neigh_avg$cum <- cumsum(TME_neigh_avg$Relabd) - TME_neigh_avg$Relabd

colorblind_dark <- c(
    "#1E1E1E",   # Blackish gray
    "#0072B2",   # Dark blue   (from the "Dark2" palette)
    "#D55E00",   # Dark orange ("Dark2" and colorblind-safe)
    "#009E73",   # Dark teal   (colorblind-friendly)
    "#332288",   # Dark purple (from Wong's palette)
    "#117733",   # Dark green  (Wong's palette)
    "#882255",   # Dark red-violet (Wong's palette)
    "#44AA99",   # Deep turquoise (Wong's palette)
    "#CC6677",   # Dusky pink (Wong's palette)
    "#AA4499"    # Dark magenta (Wong's palette)
)
relabd_barplot <- ggplot(TME_neigh_avg) +
    geom_rect(aes(xmin = cum, xmax = cum + Relabd, ymin = 0, ymax = 1, fill = Neighborhood)) +
    scale_fill_manual(values = colorblind_dark) +
    theme_void() +
    theme(legend.position="bottom")
write.csv(TME_neigh_avg, "EDA_plots/TME_avg_composition.csv")

## IHC data
IHC_data <- read.csv("data/IHC_percent.csv", row.names=1) |> na.omit() |> as.matrix()
IHC_binary <- read.csv("data/IHC_binary.csv", row.names=1) |> na.omit() |> as.matrix()


contingency_table_list <- list()
for (j in 1:3){
    for (k in (j+1):4){
        gene1 <- colnames(IHC_binary)[j]
        gene2 <- colnames(IHC_binary)[k]
        mytable <- table(IHC_binary[, j], IHC_binary[, k])
        test_result <- chisq.test(mytable)
        mytable <- as.data.frame(mytable)
        contingency_table_list[[sprintf("%s_%s", gene1, gene2)]] <- mytable
    }
}

library(openxlsx)

wb <- createWorkbook()

# Add each dataframe as a new sheet
for(name in names(contingency_table_list)) {
    addWorksheet(wb, name)
    writeData(wb, sheet = name, contingency_table_list[[name]])
}

# Save workbook to disk
saveWorkbook(wb, "EDA_plots/contingency_table.xlsx", overwrite = TRUE)

