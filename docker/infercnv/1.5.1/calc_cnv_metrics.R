args = commandArgs(trailingOnly=TRUE)

if (length(args) > 2) {
    stop("Too many arguments!")
}

filename <- args[1]
out_name <- args[2]

cnv_obj <- readRDS(filename)

# Overall CNV scores for cells.
ss_cnv <- cnv_obj@expr.data ^ 2
cnv_scores <- colSums(ss_cnv)

# Correlation with mean CNV vector within cell groups.
n_cells <- dim(cnv_obj@expr.data)[2]
corr_cnv <- rep(0, n_cells)

for (l in cnv_obj@reference_grouped_cell_indices) {
    count_mat <- as.matrix(cnv_obj@expr.data[, l])
    mean_cnv <- rowMeans(count_mat)
    corr_res <- cor(count_mat, mean_cnv)
    corr_cnv[l] <- corr_res
}

for (l in cnv_obj@observation_grouped_cell_indices) {
    count_mat <- as.matrix(cnv_obj@expr.data[, l])
    mean_cnv <- rowMeans(count_mat)
    corr_res <- cor(count_mat, mean_cnv)
    corr_cnv[l] <- corr_res
}

df <- data.frame(cbind(cnv_scores, corr_cnv))
write.table(df, file = paste0(out_name, ".csv"), quote = FALSE, sep = ",")
