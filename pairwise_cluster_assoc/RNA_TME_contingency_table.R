
rm(list=ls())

library(stringr)
library(tidyr)
library(openxlsx)



#RNA_TME
RNA_TME_result <- read.csv("multi_cluster/RNA_TME/RNA_TME_cluster.csv")
RNA_TME_result <- RNA_TME_result[, -seq(1,5)]

cnames <- colnames(RNA_TME_result)
rna_cols <- cnames[grepl("RNA", cnames)]
rna_numbers <- str_extract_all(rna_cols, "\\d+") |> unlist() |> as.integer()
TME_cols <- cnames[grepl("TME", cnames)]
TME_numbers <- str_extract_all(TME_cols, "\\d+") |> unlist() |> as.integer()

test_pvals <- expand.grid(rna_numbers, TME_numbers)
colnames(test_pvals) <- c("RNA", "TME")
test_pvals$pval <- 1
ctable <- list()
for (j in 1:nrow(test_pvals)){

    rna_name <- sprintf("RNA_%d", test_pvals$RNA[j])
    TME_name <- sprintf("TME_%d", test_pvals$TME[j])
    mytable <- table(RNA_TME_result[, rna_name], RNA_TME_result[, TME_name]) |> as.data.frame()
    mytable <- mytable %>%
        pivot_wider(
            names_from = Var2,  # Values in 'Subject' become new column names
            values_from = Freq    # Values in 'Score' populate the new columns
        )
    mytable$Var1 <- NULL
    test_result <- chisq.test(mytable,correct = TRUE)
    ctable[[sprintf("%s_%s", rna_name, TME_name)]] <- mytable
    test_pvals$pval[j] <-  test_result$p.value
}
ctable[["pvals"]] <- test_pvals
write.xlsx(
    x = ctable,
    file = "multi_cluster/RNA_TME/contingency_tables.xlsx",
    rowNames = TRUE
)
