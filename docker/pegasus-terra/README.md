# Pegasus Terra Docker Image

Pegasus is a Python package for analyzing transcriptomes of large-scale single-cell datasets. This docker image is used for interactive data analysis using Pegasus on [Terra](https://app.terra.bio) Notebook.

## Latest Release

|Image Version|Docker Image URL|Release Date|Last Modified|
|---|---|---|---|
|1.6.0|cumulusprod/pegasus-terra:1.6.0|2022/04/16|2022/04/16|

## Image Content

This image is based on terra-jupyter-base docker image. You can check out its image content [here](https://github.com/DataBiosphere/terra-docker/tree/master/terra-jupyter-base).

Besides the basic environment in terra-jupyter-base image, we also add the following main Python packages:

* [Pegasus](https://pegasus.readthedocs.io): For single-cell data analysis.
* [DemuxEM](https://demuxem.readthedocs.io): A python package on single-cell/nucleus demultiplexing.
* [Harmony-Pytorch](https://github.com/lilab-bcb/harmony-pytorch) and its dependency [PyTorch](https://pytorch.org/): Python version of Harmony algorithm on single-cell sequencing data integration using Pytorch.
* [NMF-Torch](https://github.com/lilab-bcb/nmf-torch): Pytorch implementation on Non-negative Matrix Factorization (NMF) and data integration method integrative NMF (iNMF).
* [Cirrocumulus](https://cirrocumulus.readthedocs.io): Cloud-based scalable and interactive data visualizer.
* [scVI-tools](https://scvi-tools.org/): Data integration based on Deep Learning models, and label transfer algorithms.
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html): R package on pseudo-bulk differential expression analysis.
* [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html): R package on Gene Set Enrichment Analysis (GSEA).

To see a complete list of image content, please see its [Dockerfile](https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/pegasus-terra/1.6/Dockerfile).

## Use Pegasus on Terra

To use this docker image on Terra Notebook, please refer to [this instruction](https://pegasus.readthedocs.io/en/stable/terra_notebook.html).

Tutorials on using Pegasus for interactive data analysis is [here](https://pegasus.readthedocs.io/en/stable/tutorials.html).

## Selecting Prior Versions of This Image

To use prior versions of this image, please refer to [changelog](./CHANGELOG.md) for all available versions.

## Contact Us

You can reach us by posting on our [Google Group](https://groups.google.com/forum/#!forum/cumulus-support), or [email us](mailto:cumulus-support@googlegroups.com).
