rm(list=ls())

library(stringr)
library(tidyr)
library(openxlsx)



#TME_WES
TME_WES_result <- read.csv("multi_cluster/TME_WES/TME_WES_cluster.csv")
TME_WES_result <- TME_WES_result[, -seq(1,3)]

cnames <- colnames(TME_WES_result)
TME_cols <- cnames[grepl("TME", cnames)]
TME_numbers <- str_extract_all(TME_cols, "\\d+") |> unlist() |> as.integer()
WES_cols <- cnames[grepl("WES", cnames)]
WES_numbers <- str_extract_all(WES_cols, "\\d+") |> unlist() |> as.integer()

test_pvals <- expand.grid(TME_numbers, WES_numbers)
colnames(test_pvals) <- c("TME", "WES")
test_pvals$pval <- 1
ctable <- list()
for (j in 1:nrow(test_pvals)){

    TME_name <- sprintf("TME_%d", test_pvals$TME[j])
    WES_name <- sprintf("WES_%d", test_pvals$WES[j])
    mytable <- table(TME_WES_result[, TME_name], TME_WES_result[, WES_name]) |> as.data.frame()
    mytable <- mytable %>%
        pivot_wider(
            names_from = Var2,  # Values in 'Subject' become new column names
            values_from = Freq    # Values in 'Score' populate the new columns
        )
    mytable$Var1 <- NULL
    test_result <- chisq.test(mytable,correct = TRUE)
    ctable[[sprintf("%s_%s", TME_name, WES_name)]] <- mytable
    test_pvals$pval[j] <-  test_result$p.value
}
ctable[["pvals"]] <- test_pvals
write.xlsx(
    x = ctable,
    file = "multi_cluster/TME_WES/contingency_tables.xlsx",
    rowNames = TRUE
)



