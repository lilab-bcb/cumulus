=====================================
Cumulus WDL workflows and Dockerfiles
=====================================

All of our docker images are publicly available on `Docker Hub`_ and Quay_. Our workflows use Docker Hub as the
default Docker registry. Users can use Quay as the Docker registry by entering ``quay.io/cumulus/`` for the workflow
input **"docker_registry"**, or enter a custom registry URL of their own choice.

If you use Cumulus in your research, please consider citing:

Li, B., Gould, J., Yang, Y. et al. Cumulus provides cloud-based data analysis for large-scale single-cell and
single-nucleus RNA-seq. *Nat Methods* **17**, 793–798 (2020). https://doi.org/10.1038/s41592-020-0905-x

.. _`Docker Hub`: https://cloud.docker.com/u/cumulusprod/
.. _Quay: https://quay.io/organization/cumulus

