Run SHARE-seq pipeline using shareseq_workflow
-----------------------------------------------

**shareseq_workflow** wraps 2 tools - Chromap (for ATAC-seq data) and STARsolo (for gene expression data). Chromap performs fast alignment and preprocesses chromatin profiles whereas STARsolo performs mapping and quantification for single cell RNA-seq respectively. Workflow also provides a routine to demultiplex BCL files before mapping.   

General step-by-step instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: general_steps.rst

-----------------------------------

SHARE-seq inputs and outputs 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: gex_atac.rst
