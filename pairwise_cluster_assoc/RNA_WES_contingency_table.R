
library(stringr)
library(tidyr)
library(openxlsx)

#RNA_WES
RNA_WES_result <- read.csv("multi_cluster/RNA_WES/RNA_WES_cluster.csv")
RNA_WES_result <- RNA_WES_result[, -seq(1,3)]

cnames <- colnames(RNA_WES_result)
rna_cols <- cnames[grepl("RNA", cnames)]
rna_numbers <- str_extract_all(rna_cols, "\\d+") |> unlist() |> as.integer()
wes_cols <- cnames[grepl("WES", cnames)]
wes_numbers <- str_extract_all(wes_cols, "\\d+") |> unlist() |> as.integer()

test_pvals <- expand.grid(rna_numbers, wes_numbers)
colnames(test_pvals) <- c("RNA", "WES")
test_pvals$pval <- 1
ctable <- list()
for (j in 1:nrow(test_pvals)){

    rna_name <- sprintf("RNA_%d", test_pvals$RNA[j])
    wes_name <- sprintf("WES_%d", test_pvals$WES[j])
    mytable <- table(RNA_WES_result[, rna_name], RNA_WES_result[, wes_name]) |> as.data.frame()
    mytable <- mytable %>%
        pivot_wider(
            names_from = Var2,  # Values in 'Subject' become new column names
            values_from = Freq    # Values in 'Score' populate the new columns
        )
    mytable$Var1 <- NULL
    test_result <- chisq.test(mytable,correct = TRUE)
    ctable[[sprintf("%s_%s", rna_name, wes_name)]] <- mytable
    test_pvals$pval[j] <-  test_result$p.value
}
ctable[["pvals"]] <- test_pvals
write.xlsx(
    x = ctable,
    file = "multi_cluster/RNA_WES/contingency_tables.xlsx",
    rowNames = TRUE
)





