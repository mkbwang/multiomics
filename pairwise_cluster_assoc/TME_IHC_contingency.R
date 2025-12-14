
rm(list=ls())
library(stringr)
library(tidyr)
library(openxlsx)

#TME_IHC
TME_IHC_result <- read.csv("multi_cluster/TME_IHC/TME_IHC_cluster.csv")
TME_IHC_result <- TME_IHC_result[, -seq(1, 3)]

cnames <- colnames(TME_IHC_result)
TME_cols <- cnames[grepl("TME", cnames)]
TME_numbers <- str_extract_all(TME_cols, "\\d+") |> unlist() |> as.integer()
IHC_cols <- cnames[-seq(1,3)]

test_pvals <- expand.grid(TME_numbers, IHC_cols)
colnames(test_pvals) <- c("TME", "IHC")
test_pvals$pval <- 1
ctable <- list()
for (j in 1:nrow(test_pvals)){

    TME_name <- sprintf("TME_%d", test_pvals$TME[j])
    IHC_name <- test_pvals$IHC[j] |> as.character()
    mytable <- table(TME_IHC_result[, TME_name], TME_IHC_result[, IHC_name]) |> as.data.frame()
    mytable <- mytable %>%
        pivot_wider(
            names_from = Var2,  # Values in 'Subject' become new column names
            values_from = Freq    # Values in 'Score' populate the new columns
        )
    mytable$Var1 <- NULL
    test_result <- chisq.test(mytable,correct = TRUE)
    ctable[[sprintf("%s_%s", TME_name, IHC_name)]] <- mytable
    test_pvals$pval[j] <-  test_result$p.value
}
ctable[["pvals"]] <- test_pvals
write.xlsx(
    x = ctable,
    file = "multi_cluster/TME_IHC/contingency_tables.xlsx",
    rowNames = TRUE
)







