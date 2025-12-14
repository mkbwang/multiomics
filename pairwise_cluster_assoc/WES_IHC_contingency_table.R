
rm(list=ls())
library(stringr)
library(tidyr)
library(openxlsx)

#WES_IHC
WES_IHC_result <- read.csv("multi_cluster/WES_IHC/WES_IHC_cluster.csv")
WES_IHC_result <- WES_IHC_result[, -1]

cnames <- colnames(WES_IHC_result)
WES_cols <- cnames[grepl("WES", cnames)]
WES_numbers <- str_extract_all(WES_cols, "\\d+") |> unlist() |> as.integer()
IHC_cols <- cnames[-seq(1,4)]

test_pvals <- expand.grid(WES_numbers, IHC_cols)
colnames(test_pvals) <- c("WES", "IHC")
test_pvals$pval <- 1
ctable <- list()
for (j in 1:nrow(test_pvals)){

    WES_name <- sprintf("WES_%d", test_pvals$WES[j])
    IHC_name <- test_pvals$IHC[j] |> as.character()
    mytable <- table(WES_IHC_result[, WES_name], WES_IHC_result[, IHC_name]) |> as.data.frame()
    mytable <- mytable %>%
        pivot_wider(
            names_from = Var2,  # Values in 'Subject' become new column names
            values_from = Freq    # Values in 'Score' populate the new columns
        )
    mytable$Var1 <- NULL
    test_result <- chisq.test(mytable,correct = TRUE)
    ctable[[sprintf("%s_%s", WES_name, IHC_name)]] <- mytable
    test_pvals$pval[j] <-  test_result$p.value
}
ctable[["pvals"]] <- test_pvals
write.xlsx(
    x = ctable,
    file = "multi_cluster/WES_IHC/contingency_tables.xlsx",
    rowNames = TRUE
)







