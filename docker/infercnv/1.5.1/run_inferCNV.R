library(optparse)
library(infercnv)

option.list <- list(
    make_option("--sample-id",
                action = "store",
                type = "character",
                default = NULL,
                dest = "sample_id",
                metavar = "Sample Name",
                help = paste("Sample Name.",
                             "[Default %default]")),
    make_option("--raw-counts-matrix",
                action = "store",
                type = "character",
                default = NULL,
                dest = "raw_counts_matrix",
                metavar = "Raw Counts Expression Data",
                help = paste("the matrix of genes (rows) vs. cells (columns) ",
                             "containing the raw counts. It'll be read via read.table()",
                             "[Default %default]")),
    make_option("--gene-order-file",
                action = "store",
                type = "character",
                default = NULL,
                dest = "gene_order_file",
                metavar = "Gene Order File",
                help = paste("data file containing the positions of ",
                             "each gene along each chromosome in the genome. ",
                             "[Default %default]")),
    make_option("--annotations-file",
                action = "store",
                type = "character",
                default = NULL,
                dest = "annotations_file",
                metavar = "Annotation File",
                help = paste("a description of the cells, indicating ",
                             "the cell type classifications. ",
                             "[Default %default]")),
    make_option("--protocol",
                action = "store",
                type = "character",
                default = "tenX",
                dest = "protocol",
                metavar = "Protocol",
                help = paste("Cutoff is chosen depending on Protocols.",
                             "Available options are 'tenX' and 'SMART-Seq2'",
                             "[Default %default]")),
    make_option("--ref-group-names",
                action = "store",
                type = "character",
                default = NULL,
                dest = "ref_group_names",
                metavar = "Reference Groups Names",
                help = paste("Names of groups from raw_counts_matrix whose cells",
                            "are to be used as reference groups.",
                            "[Default %default]")),
    make_option("--out-dir",
                action = "store",
                type = "character",
                default = ".",
                dest = "out_dir",
                metavar = "Output Directory",
                help = paste("Path to directory to deposit outputs. ",
                            "[Default %default]")),
    make_option("--cluster-by-groups",
                action = "store_true",
                type = "logical",
                default = FALSE,
                dest = "cluster_by_groups",
                metavar = "Cluster by Groups",
                help = paste("If observations are defined according to groups ",
                            "(ie. patients), each group of cells will be ",
                            "clustered separately. ([Default %default]",
                            ", instead will use k_obs_groups setting)")),
    make_option("--HMM",
                action = "store_true",
                type = "logical",
                default = FALSE,
                dest = "HMM",
                metavar = "HMM",
                help = paste("when set to True, runs HMM to predict CNV level. ",
                             "[Default %default]"))
)

args <- parse_args(OptionParser(option_list = option.list))

if (args$protocol == "tenX") {
    cutoff <- 0.1
} else {
    cutoff <- 1
}

if (!is.null(args$ref_group_names)) {
    args$ref_group_names <- strsplit(args$ref_group_names, ",")[[1]]
} else {
    args$ref_group_names <- c()
}

cnv_obj <- CreateInfercnvObject(raw_counts_matrix = args$raw_counts_matrix,
                                annotations_file = args$annotations_file,
                                delim = '\t',
                                gene_order_file = args$gene_order_file,
                                ref_group_names = args$ref_group_names,
                                min_max_counts_per_cell = c(-Inf, +Inf))

cnv_obj <- infercnv::run(cnv_obj, cutoff = cutoff, out_dir = args$out_dir, cluster_by_groups = args$cluster_by_groups, HMM = args$HMM)
