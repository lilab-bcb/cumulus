if (!'loomR' %in% rownames(installed.packages())) {
    require("devtools")
    devtools::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
}

library(loomR)
library(Seurat)

convert_loom_to_seurat <- function(from) {
    isV3 <- startsWith(toString(packageVersion("Seurat")), "3.")

    if (!isV3) {
        stop("Convertion from loom to Seurat object only works for Seurat v3!")
    }

    x.loom <- connect(from, mode = 'r')
    seurat.object <- as.Seurat(x.loom, cells = 'obs_names', features = 'var_names')

    cell.attrs <- hdf5r::list.datasets(
        object = x.loom,
        path = 'col_attrs',
        full.names = FALSE,
        recursive = FALSE
    )

    embed.attrs <- Filter(
        f = function(embed) {
            return(length(x.loom[["col_attrs"]][[embed]]$dims) == 2)
        },
        x = cell.attrs
    )

    # For cell embeddings.
    for (embed.key in embed.attrs) {
        if (startsWith(embed.key, "X_")) {
            embed.type <- substr(embed.key, 3, nchar(embed.key))
            embed.name <- switch(embed.type, "pca" = "PC", "tsne" = "tSNE", "diffmap_pca" = "DIFFMAPPCA", toupper(embed.type))
            embed.name <- paste0(embed.name, "_")
            embed.content <- t(x = x.loom[['col_attrs']][[embed.key]][, ])
            colnames(x = embed.content) <- paste0(embed.name, 1:ncol(x = embed.content))
            rownames(x = embed.content) <- colnames(x = seurat.object)
            dim.reduc <- CreateDimReducObject(embeddings = embed.content, assay = "RNA", key = embed.name)
            seurat.object[[embed.type]] <- dim.reduc
            
            message(paste0(embed.key, " is added."))
        }
    }
    return(seurat.object)
}
