To utilize cellranger arc/cellranger multi/cellranger count for single-cell multiomics, follow the specific instructions below. In particular, we put each single modality in one separate lin in the sample sheet as described above. We then use the *Link* column to link multiple modalities together. Depending on the modalities included, *cellranger arc* (Multiome ATAC + Gene Expression), *cellranger multi* (CellPlex), or *cellranger count* (Feature Barcode) will be triggered. Note that cumulus_feature_barcoding/demuxEM would not be triggered for hashing/citeseq in this setting.

Sample sheet
++++++++++++

#. **Reference** column.

	Pre-built Multiome ATAC + Gene Expression references are summarized below. CellPlex and Feature Barcode use the same reference as in Single-cell and single-nucleus RNA-seq.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCh38-2020-A_arc_v2.0.0**
		  - Human GRCh38 sequences (GENCODE v32/Ensembl 98), cellranger arc reference 2.0.0
		* - **mm10-2020-A_arc_v2.0.0**
		  - Mouse GRCm38 sequences (GENCODE vM23/Ensembl 98), cellranger arc reference 2.0.0

#. *DataType* column.

	For each modality, set it to the corresponding data type.

#. *FetureBarcodeFile* column.

	For RNA-seq modality, only set this if a target panel is provided.
	For CMO (CellPlex), provide sample name - CMO tag association as follows::

		sample1,CMO301|CMO302
		sample2,CMO303

	For CITESeq, Perturb-seq and hashing, provide one CSV file as defined in `Feature Barcode Reference`_. Note that one feature barcode reference should be provided for all feature-barcode related modalities (e.g. *citeseq*, *hashing*, *crispr*) and all these modalities should put the same reference file in *FeatureBarcodeFile* column.

#. *Link* column.

	Put a sample unique link name for all modalities that are linked.

#. Example::

	Sample,Reference,Flowcell,DataType,FeatureBarcodeFile,Link
	sample1_rna,GRCh38-2020-A_arc_v2.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ/Fastq,rna,,sample1
	sample1_atac,GRCh38-2020-A_arc_v2.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ/Fastq,atac,,sample1
	sample2_rna,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZX/Fastq,,rna,,sample2
	sample2_cmo,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZX/Fastq,cmo,gs://fc-e0000000-0000-0000-0000-000000000000/cmo.csv,sample2
	sample3_rna,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZY/Fastq,rna,,sample3
	sample3_citeseq,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZY/Fastq,citeseq,gs://fc-e0000000-0000-0000-0000-000000000000/feature_ref.csv,sample3


In the above example, three linked samples are provided. *cellranger arc*, *cellranger multi* and *cellranger count* will be triggered respectively.

Workflow input
++++++++++++++

For single-cell multiomics data, ``cellranger_workflow`` takes sequencing reads as input (FASTQ files, or TAR files containing FASTQ files), and runs ``cellranger-arc ount/cellranger multi/cellranger count``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

.. list-table::
	:widths: 5 30 30 20
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_csv_file**
	  - Sample Sheet (contains Sample, Reference, Flowcell as required and Chemistry, DataType, FeatureBarcodeFile, Link as optional)
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
	  -
	* - **output_directory**
	  - Output directory
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
	  -
	* - include_introns
	  - Turn this option on to also count reads mapping to intronic regions. With this option, users do not need to use pre-mRNA references. Note that if this option is set, cellranger_version must be >= 5.0.0. This option is used by *cellranger multi* and *cellranger count*.
	  - true
	  - true
	* - no_bam
	  - Turn this option on to disable BAM file generation. This option is only available if cellranger_version >= 5.0.0. This option is used by *cellranger-arc count*, *cellranger multi* and *cellranger count*.
	  - false
	  - false
	* - expect_cells
	  - Expected number of recovered cells. Mutually exclusive with force_cells. This option is used by *cellranger multi* and *cellranger count*.
	  - 3000
	  -
	* - force_cells
	  - Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells. This option is used by *cellranger multi* and *cellranger count*.
	  - 6000
	  -
	* - arc_gex_exclude_introns
	  - | Disable counting of intronic reads. In this mode, only reads that are exonic and compatible with annotated splice junctions in the reference are counted.
	    | **Note:** using this mode will reduce the UMI counts in the feature-barcode matrix.
	  - false
	  - false
	* - arc_min_atac_count
	  - | Cell caller override to define the minimum number of ATAC transposition events in peaks (ATAC counts) for a cell barcode.
	    | **Note:** this input must be specified in conjunction with ``arc_min_gex_count`` input.
	    | With both inputs set, a barcode is defined as a cell if it contains at least ``arc_min_atac_count`` ATAC counts AND at least ``arc_min_gex_count`` GEX UMI counts.
	  - 100
	  -
	* - arc_min_gex_count
	  - | Cell caller override to define the minimum number of GEX UMI counts for a cell barcode.
	    | **Note:** this input must be specified in conjunction with ``arc_min_atac_count``. See the description of ``arc_min_atac_count`` input for details.
	  - 200
	  -
	* - peaks
	  - A 3-column BED file of peaks to override cellranger arc peak caller. Peaks must be sorted by position and not contain overlapping peaks; comment lines beginning with ``#`` are allowed
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/common_peaks.bed"
	  -
	* - secondary
	  - Perform Cell Ranger secondary analysis (dimensionality reduction, clustering, etc.). This option is used by *cellranger multi* and *cellranger count*.
	  - false
	  - false
	* - cmo_set
	  - CMO set CSV file, delaring CMO constructs and associated barcodes. See `CMO reference`_ for details. Used only for *cellranger multi*.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/cmo_set.csv"
	  -
	* - cellranger_version
	  - cellranger version, could be: 9.0.1, 9.0.0, 8.0.1, 8.0.0, 7.2.0, 7.1.0, 7.0.1, 7.0.0
	  - "9.0.1"
	  - "9.0.1"
	* - cellranger_arc_version
	  - cellranger-arc version, could be: ``2.0.2.strato`` (compatible with workflow v2.6.1+), ``2.0.2.custom-max-cell`` (with max_cell threshold set to 80,000), ``2.0.2`` (compatible with workflow v2.6.0 or earlier), ``2.0.1``, ``2.0.0``
	  - "2.0.2.strato"
	  - "2.0.2.strato"
	* - docker_registry
	  - Docker registry to use for cellranger_workflow. Options:

	  	- "quay.io/cumulus" for images on Red Hat registry;

	  	- "cumulusprod" for backup images on Docker Hub.
	  - "quay.io/cumulus"
	  - "quay.io/cumulus"
	* - acronym_file
	  - | The link/path of an index file in TSV format for fetching preset genome references, chemistry barcode inclusion lists, etc. by their names.
	    | Set an GS URI if *backend* is ``gcp``; an S3 URI for ``aws`` backend; an absolute file path for ``local`` backend.
	  - "s3://xxxx/index.tsv"
	  - "gs://cumulus-ref/resources/cellranger/index.tsv"
	* - zones
	  - Google cloud zones. For GCP Batch backend, the zones are automatically restricted by the Batch settings.
	  - "us-central1-a us-west1-a"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - num_cpu
	  - Number of cpus to request for one node for cellranger vdj
	  - 32
	  - 32
	* - memory
	  - Memory size string for cellranger and cellranger vdj
	  - "120G"
	  - "120G"
	* - count_disk_space
	  - Disk space in GB needed for cellranger count
	  - 500
	  - 500
	* - arc_num_cpu
	  - Number of cpus to request for one node for cellranger-arc count
	  - 64
	  - 64
	* - arc_memory
	  - Memory size string for cellranger-arc count
	  - "160G"
	  - "160G"
	* - arc_disk_space
	  - Disk space in GB needed for cellranger-arc count
	  - 700
	  - 700
	* - backend
	  - Cloud backend for file transfer. Available options:

	  	- "gcp" for Google Cloud;
	  	- "aws" for Amazon AWS;
	  	- "local" for local machine.
	  - "gcp"
	  - "gcp"
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2
	* - awsQueueArn
	  - The AWS ARN string of the job queue to be used. This only works for ``aws`` backend.
	  - "arn:aws:batch:us-east-1:xxx:job-queue/priority-gwf"
	  - ""

Workflow output
+++++++++++++++

See the table below for important sc/snRNA-seq outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - cellranger_arc_count.output_count_directory / cellranger_multi.output_multi_directory / cellranger_count_fbc.output_count_directory
	  - Array[String]
	  - Subworkflow output. A list of cloud urls containing *cellranger-arc count*, *cellranger multi* or *cellranger count* outputs, one url per sample.
	* - cellranger_arc_count.output_web_summary / cellranger_count_fbc.output_web_summary
	  - Array[File]
	  - A list of htmls visualizing QCs for each sample (*cellranger-arc count* / *cellranger count* output).
	* - collect_summaries_arc.metrics_summaries / collect_summaries_fbc.metrics_summaries
	  - File
	  - A excel spreadsheet containing QCs for each sample.


.. _Feature Barcode Reference: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#feature-ref
.. _CMO reference: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#cmoreference
