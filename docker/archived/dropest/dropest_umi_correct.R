#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
library(dropestr)
library(parallel)
count_matrix <- readRDS(args[1])
reads_per_umi_per_cell <- count_matrix$reads_per_umi_per_cell
cores <- as.integer(args[2]) #detectCores(logical = FALSE)
count_matrix <- CorrectUmiSequenceErrors(reads_per_umi_per_cell, mc.cores = cores, verbosity.level = 2)
saveRDS(count_matrix, "umi_corrected.rds")


# umis_distribution <- GetUmisDistribution(reads_per_umi_per_cell$reads_per_umi)
# umi_probabilities <- umis_distribution / sum(umis_distribution)
# PlotUmisDistribution(reads_per_umi_per_cell$reads_per_umi)
#
#
# umis_per_gene <- sapply(reads_per_umi_per_cell$reads_per_umi, length)
# max_umi_per_gene <- max(umis_per_gene)
# collision_info <- FillCollisionsAdjustmentInfo(umi_probabilities, max_umi_per_gene, verbose = T)
# c(max_umi_per_gene, collision_info[max_umi_per_gene])
#
# correction_info <- PrepareUmiCorrectionInfo(umi_probabilities, collision_info[max_umi_per_gene],
# quants.num = 50, verbosity.level = 2)
#
#
# cm <- CorrectUmiSequenceErrors(reads_per_umi_per_cell, umi.probabilities = umi_probabilities, collisions.info = collision_info,
# correction.info = correction_info, mc.cores = 5, verbosity.level = 2)
# saveRDS(cm, "umi_corrected.rds")
