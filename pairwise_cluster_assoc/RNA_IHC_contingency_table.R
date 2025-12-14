
rm(list=ls())
library(stringr)
library(tidyr)
library(openxlsx)

#RNA_IHC
RNA_IHC_result <- read.csv("multi_cluster/RNA_IHC/RNA_IHC_cluster.csv")
RNA_IHC_result <- RNA_IHC_result[, -seq(1,3)]

cnames <- colnames(RNA_IHC_result)
rna_cols <- cnames[grepl("RNA", cnames)]
rna_numbers <- str_extract_all(rna_cols, "\\d+") |> unlist() |> as.integer()
IHC_cols <- cnames[-seq(1,3)]

test_pvals <- expand.grid(rna_numbers, IHC_cols)
colnames(test_pvals) <- c("RNA", "IHC")
test_pvals$pval <- 1
ctable <- list()
for (j in 1:nrow(test_pvals)){

    rna_name <- sprintf("RNA_%d", test_pvals$RNA[j])
    IHC_name <- test_pvals$IHC[j] |> as.character()
    mytable <- table(RNA_IHC_result[, rna_name], RNA_IHC_result[, IHC_name]) |> as.data.frame()
    mytable <- mytable %>%
        pivot_wider(
            names_from = Var2,  # Values in 'Subject' become new column names
            values_from = Freq    # Values in 'Score' populate the new columns
        )
    mytable$Var1 <- NULL
    test_result <- chisq.test(mytable,correct = TRUE)
    ctable[[sprintf("%s_%s", rna_name, IHC_name)]] <- mytable
    test_pvals$pval[j] <-  test_result$p.value
}
ctable[["pvals"]] <- test_pvals
write.xlsx(
    x = ctable,
    file = "multi_cluster/RNA_IHC/contingency_tables.xlsx",
    rowNames = TRUE
)







