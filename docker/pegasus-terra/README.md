# Pegasus Terra Docker Image

Pegasus is a Python package for analyzing transcriptomes of large-scale single-cell datasets. This docker image is used for interactive data analysis using Pegasus on [Terra](https://app.terra.bio) Notebook.

## Latest Release

|Image Version|Docker Image URL|Release Date|
|---|---|---|
|0.16.11|cumulusprod/pegasus-terra:0.16.11|2020/03/02|

## Image Content

This image is based on terra-jupyter-base docker image. You can check out its image content [here](https://github.com/DataBiosphere/terra-docker/tree/master/terra-jupyter-base).

Besides the basic environment in terra-jupyter-base image, we also add the following main Python packages:

* [Pegasus](https://pegasus.readthedocs.io): for single-cell data analysis.
* [scPlot](https://github.com/klarman-cell-observatory/scPlot): Interactive single-cell plotting library used by Pegasus.
* [Harmony-Pytorch](https://github.com/lilab-bcb/harmony-pytorch) and its dependency [PyTorch](https://pytorch.org/): Python version of Harmony algorithm on single-cell sequencing data integration.

To see a complete list of image content, please see its [Dockerfile](https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/pegasus-terra/0.16.8/Dockerfile).

## Use Pegasus on Terra

To use this docker image on Terra Notebook, please refer to [this instruction](https://pegasus.readthedocs.io/en/latest/terra_notebook.html).

Tutorials on using Pegasus for interactive data analysis is [here](https://pegasus.readthedocs.io/en/latest/tutorials.html).

## Selecting Prior Versions of This Image

To use prior versions of this image, please refer to [changelog](./CHANGELOG.md) for all available versions.

## Contact Us

You can reach us by posting on our [Google Group](https://groups.google.com/forum/#!forum/cumulus-support), or [email us](mailto:cumulus-support@googlegroups.com).