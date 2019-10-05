### This script is modified from conversion.R of Seurat V2 source code. Seurat is licensed under GPL 3.0

#' @describeIn Convert from Anndata file to a Seurat object
#' @importFrom reticulate py_to_r
#' @export Seurat, Matrix, reticulate
#' @method convert_h5ad_to_seurat
#'
#'

library(Seurat)
library(Matrix)
library(reticulate)

convert_h5ad_to_seurat <- function(from) {
	isV3 <- startsWith(toString(packageVersion("Seurat")), "3.")

	meta.data <- py_to_r(from$obs)
	for (key in colnames(meta.data)) {
		if (from$obs[key]$dtype$name == "category") {
			meta.data[key] = py_to_r(from$obs[key]$astype("str"))
		}
	}
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

	meta.features <- py_to_r(from$var)
	for (key in colnames(meta.features)) {
		if (from$var[key]$dtype$name == "category") {
			meta.features[key] = py_to_r(from$var[key]$astype("str"))
		}
	}

	raw.data.matrix <- sparseMatrix(
		i = as.numeric(x = from$raw$X$indices),
		p = as.numeric(x = from$raw$X$indptr),
		x = as.numeric(x = from$raw$X$data),
		index1 = FALSE
	)
	rownames(x = raw.data.matrix) <- rownames(x = py_to_r(from$raw$var))
	colnames(x = raw.data.matrix) <- rownames(x = meta.data)

	if (isV3) {
		seurat.object <- CreateSeuratObject(
			counts = raw.data.matrix,
			meta.data = meta.data,
		)
		tmp = seurat.object[["RNA"]]
		tmp[[names(x = meta.features)]] = meta.features
		seurat.object[["RNA"]] = tmp
	} else {
		seurat.object <- CreateSeuratObject(
			raw.data = raw.data.matrix,
			meta.data = meta.data
		)
	}

	data.matrix <- sparseMatrix(
		i = as.numeric(x = from$X$indices),
		p = as.numeric(x = from$X$indptr),
		x = as.numeric(x = from$X$data),
		index1 = FALSE
	)
	rownames(x = data.matrix) <- rownames(x = meta.features)
	colnames(x = data.matrix) <- rownames(x = meta.data)

	if (isV3) {
		seurat.object <- SetAssayData(
			object = seurat.object,
			slot = "data",
			new.data = data.matrix
		)
	} else {
		seurat.object <- SetAssayData(
			object = seurat.object,
			assay.type = "RNA",
			slot = "data",
			new.data = data.matrix
		)
	}

	scale.data.matrix <- t(x = py_to_r(from$uns["scale.data"]))
	rownames(x = scale.data.matrix) <- py_to_r(from$uns["scale.data.rownames"])
	colnames(x = scale.data.matrix) <- rownames(x = meta.data)

	if (isV3) {
		seurat.object <- SetAssayData(
			object = seurat.object,
			slot = "scale.data",
			new.data = scale.data.matrix
		)
	} else {
		seurat.object <- SetAssayData(
			object = seurat.object,
			assay.type = "RNA",
			slot = "scale.data",
			new.data = scale.data.matrix
		)
	}

	obsm_keys <- toString(from$obsm$keys())
	obsm_keys <- gsub("KeysView(AxisArrays with keys: ", "", obsm_keys, fixed = TRUE)
	obsm_keys <- substr(obsm_keys, 1, nchar(obsm_keys) - 1)
	obsm_keys <- strsplit(obsm_keys, split = ", ", fixed = TRUE)[[1]]

	for (embed.key in obsm_keys) {
		if (startsWith(embed.key, "X_")) {
			embed.type <- substr(embed.key, 3, nchar(embed.key))
			embed.name <- switch(embed.type, "pca" = "PC", "tsne" = "tSNE", toupper(embed.type))
			embed.name <- paste0(embed.name, "_")
			if (!isV3 && embed.name == "PC_") { embed.name = "PC" }
			embed.content = py_to_r(from$obsm$get(embed.key))
			colnames(x = embed.content) <- paste0(embed.name, 1:ncol(x = embed.content))
			rownames(x = embed.content) <- rownames(x = meta.data)
			if (isV3) {
				dim.reduc <- CreateDimReducObject(embeddings = embed.content, assay = "RNA", key = embed.name)
				seurat.object[[embed.type]] <- dim.reduc
			} else {
				seurat.object <- SetDimReduction(
					object = seurat.object,
					reduction.type = embed.type,
					slot = "cell.embeddings",
					new.data = embed.content
				)
				seurat.object <- SetDimReduction(
					object = seurat.object,
					reduction.type = embed.type,
					slot = "key",
					new.data = embed.name
				)
			}
		}
	}

	return(seurat.object)
}
