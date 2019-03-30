#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
library(dropestr)
library(parallel)
count_matrix <- readRDS(args[1])
cores <- detectCores(logical = FALSE)
count_matrix <- CorrectUmiSequenceErrors(count_matrix$reads_per_umi_per_cell, mc.cores = cores)
saveRDS(count_matrix, "umi_corrected.rds")
