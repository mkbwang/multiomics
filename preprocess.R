
library(readxl)
# read RNAseq file
## raw counts
rna_counts <- read_excel("excel_data/RNA_normalized_protein_coding.xlsx", sheet=1) |>
    as.data.frame()
rownames(rna_counts) <- rna_counts[, 1]
rna_counts[, 1] <- NULL
rna_counts <- round(rna_counts)
## do my own normalization with DESEQ2
library(DESeq2)
colData <- data.frame(row.names = colnames(rna_counts),
                      group = factor(rep("dummy", ncol(rna_counts))))
dds <- DESeqDataSetFromMatrix(countData=rna_counts, colData=colData, design = ~1)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)
norm_counts <- assay(vsd)

rna_normalized <- read_excel("excel_data/RNA_normalized_protein_coding.xlsx", sheet=3) |>
    as.data.frame()
rownames(rna_normalized) <- rna_normalized[, 1]
rna_normalized[, 1] <- NULL

rna_counts <- rna_counts[rownames(rna_normalized), ]
norm_counts <- norm_counts[rownames(rna_normalized), ]
all(abs(rna_normalized - norm_counts) < 1e-10) # confirm that the normalization comes from DESeq2
write.csv(rna_normalized, "RNA_normalized.csv", quote=FALSE)
write.csv(rna_counts, "data/RNA_counts.csv", quote=FALSE)

# read whole exome sequencing
wes_data <- read_excel("excel_data/WES_filtered.VAF5.binary.xlsx") |> as.data.frame()
rownames(wes_data) <- wes_data[, 1]
wes_data[, 1] <- NULL
write.csv(wes_data, "WES_data.csv", quote=FALSE)

# read IHC data
IHC_data <- read_excel("excel_data/IHC_Sept_2025.xlsx") |> as.data.frame()
rownames(IHC_data) <- IHC_data[, 1]
IHC_data[, 1] <- NULL
feature_names <- colnames(IHC_data)
IHC_data_score <- IHC_data[, c("T-bet score", "GATA-3 score", "CXCR3 score",
                               "CCR4 score")]
colnames(IHC_data_score) <- c("T-bet", "GATA-3", "CXCR3",
                              "CCR4")
write.csv(IHC_data_score, "IHC_binary.csv", quote=FALSE)


IHC_data_percent <- IHC_data[, c("T-bet %postive", "GATA3 %positive",
                                 "CXCR3 % pos", "CCR4%")]
colnames(IHC_data_percent) <- c("T-bet", "GATA-3", "CXCR3", "CCR4")
write.csv(IHC_data_percent, "IHC_percent.csv", quote=FALSE)


# read TME data

TME <- read_excel("excel_data/TME.xlsx") |> as.data.frame()
rownames(TME) <- TME[, 1]
TME[, 1] <- NULL
TME <- TME[, seq(1, 12)] |> na.omit()
colnames(TME) <- c(sprintf("Neigh%d", seq(0, 9)), "Cytotoxic", "Macrophage")
write.csv(TME, "TME.csv", quote=FALSE)


rna_samples <- colnames(rna_normalized)
wes_samples <- colnames(wes_data)
IHC_samples <- rownames(IHC_data_score)

length(intersect(IHC_samples,wes_samples))
