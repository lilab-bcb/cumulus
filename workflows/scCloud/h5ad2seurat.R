### This script is modified from conversion.R of Seurat's github. Seurat is licensed under GPL 3.0


#' @param X.slot Seurat slot to transfer anndata X into. Default is scale.data
#' @param raw.slot Seurat slot to transfer anndata raw into. Default is data
#' @describeIn Convert from Anndata file to a Seurat object
#' @importFrom reticulate py_to_r
#' @export
#' @method Convert anndata.base.AnnData
#'
#'
Convert.anndata.base.AnnData <- function(
	from,
	to,
	...
) {
	object.to <- switch(
		EXPR = to,
		'seurat' = {
			raw.data.matrix <- sparseMatrix(
				i = as.numeric(x = from$raw$X$indices),
				p = as.numeric(x = from$raw$X$indptr),
				x = as.numeric(x = from$raw$X$data),
				index1 = FALSE
			)
			rownames(x = raw.data.matrix) <- rownames(x = py_to_r(from$var))
			colnames(x = raw.data.matrix) <- rownames(x = py_to_r(from$obs))
			meta.data <- py_to_r(from$obs)
			if ("n_counts" %in% colnames(x = meta.data)) {
				colnames(x = meta.data) <- gsub(
					pattern = "n_counts",
					replacement = "nUMI",
					x = colnames(x = meta.data)
				)
			}
			if ("n_genes" %in% colnames(x = meta.data)) {
				colnames(x = meta.data) <- gsub(
					pattern = "n_genes",
					replacement = "nGene",
					x = colnames(x = meta.data)
				)
			}
			seurat.object <- CreateSeuratObject(
				raw.data = raw.data.matrix,
				meta.data = meta.data
			)

			data.matrix <- sparseMatrix(
				i = as.numeric(x = from$X$indices),
				p = as.numeric(x = from$X$indptr),
				x = as.numeric(x = from$X$data),
				index1 = FALSE
			)
			rownames(x = data.matrix) <- rownames(x = py_to_r(from$var))
			colnames(x = data.matrix) <- rownames(x = py_to_r(from$obs))
			seurat.object <- SetAssayData(
				object = seurat.object,
				assay.type = "RNA",
				slot = "data",
				new.data = data.matrix
			)

			scale.data.matrix <- t(x = py_to_r(from$uns["scale.data"]))
			rownames(x = scale.data.matrix) <- py_to_r(from$uns["scale.data.rownames"])
			colnames(x = scale.data.matrix) <- rownames(x = py_to_r(from$obs))
			seurat.object <- SetAssayData(
				object = seurat.object,
				assay.type = "RNA",
				slot = "scale.data",
				new.data = scale.data.matrix
			)

			dr.dict <- list(tSNE_ = "tsne", PC = "pca", UMAP = "umap", FLE = "fle", DIFFMAP_PCA = "diffmap_pca")
			drs <- unlist(x = py_to_r(from$obsm$keys()))
			for (dr.key in names(dr.dict)) {
				dr.name <- dr.dict[[dr.key]]
				dr.pkey <- paste0("X_", dr.name)
				if (dr.pkey %in% drs) {
					dr.embed <- py_to_r(from$obsm[[dr.pkey]])
					colnames(x = dr.embed) <- paste0(dr.key, 1:ncol(x = dr.embed))
					rownames(x = dr.embed) <- seurat.object@cell.names
					seurat.object <- SetDimReduction(
						object = seurat.object,
						reduction.type = dr.name,
						slot = "cell.embeddings",
						new.data = dr.embed
					)
					seurat.object <- SetDimReduction(
						object = seurat.object,
						reduction.type = dr.name,
						slot = "key",
						new.data = dr.key
					)
				}
			}
			seurat.object
		},
		stop(paste0("Cannot convert AnnData objects to class '", to, "'"))
	)
	return(object.to)
}
