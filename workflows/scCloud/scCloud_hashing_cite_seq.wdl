import "https://api.firecloud.org/ga4gh/v1/tools/scCloud:tasks/versions/23/plain-WDL/descriptor" as tasks
# import "../scCloud/scCloud_tasks.wdl" as tasks

workflow scCloud_hashing_cite_seq {
	# An sample sheet contains RNA and ADT data correspondence
	File input_sample_sheet
	# Output directory, gs url
	String output_directory

	# scCloud version, default to "0.8.0"
	String? sccloud_version = "0.8.0"
	# Google cloud zones, default to "us-east1-d us-west1-a us-west1-b"
	String? zones = "us-east1-d us-west1-a us-west1-b"
	# Number of cpus
	Int? num_cpu = 8
	# Memory string
	String? memory = "10G"
	# Disk space in GB
	Int? disk_space = 20
	# Number of preemptible tries 
	Int? preemptible = 2

	# Reference genome
	String? genome

	# demuxEM parameters
	# Only demultiplex cells/nuclei with at least <demuxEM_min_num_genes> expressed genes. [default: 100]
	Int? demuxEM_min_num_genes
	# The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse.
	Float? demuxEM_alpha_on_samples
	# Only demultiplex cells/nuclei with at least <demuxEM_min_num_umis> of UMIs. [default: 100]
	Int? demuxEM_min_num_umis
	# Any cell/nucleus with less than <count> hashtags from the signal will be marked as unknown. [default: 10.0]
	Float? demuxEM_min_signal_hashtag
	# The random seed used in the KMeans algorithm to separate empty ADT droplets from others. [default: 0]
	Int? demuxEM_random_state
	# Generate a series of diagnostic plots, including the background/signal between HTO counts, estimated background probabilities, HTO distributions of cells and non-cells etc. [default: true]
	Boolean? demuxEM_generate_diagnostic_plots = true
	# Generate violin plots using gender-specific genes (e.g. Xist). <demuxEM_generate_gender_plot> is a comma-separated list of gene names.
	String? demuxEM_generate_gender_plot

	# merge_rna_adt parameters
	# A CSV file containing the IgG control information for each antibody.
	File? antibody_control_csv



	call tasks.generate_hashing_cite_seq_tasks as generate_hashing_cite_seq_tasks {
		input:
			input_sample_sheet = input_sample_sheet,
			sccloud_version = sccloud_version,
			zones = zones,
			preemptible = preemptible
	}

	if (generate_hashing_cite_seq_tasks.hashing_ids[0] != '') {
		scatter (hashing_id in generate_hashing_cite_seq_tasks.hashing_ids) {
			call tasks.run_scCloud_demuxEM as run_scCloud_demuxEM {
				input:
					input_adt_csv = generate_hashing_cite_seq_tasks.id2adt[hashing_id],
					input_raw_gene_bc_matrices_h5 = generate_hashing_cite_seq_tasks.id2rna[hashing_id],
					output_dir = sub(output_directory, "/+$", ""),
					output_name = hashing_id,
					genome = genome,
					alpha_on_samples = demuxEM_alpha_on_samples,
					min_num_genes = demuxEM_min_num_genes,
					min_num_umis = demuxEM_min_num_umis,
					min_signal_hashtag = demuxEM_min_signal_hashtag,
					random_state = demuxEM_random_state,
					generate_diagnostic_plots = demuxEM_generate_diagnostic_plots,
					generate_gender_plot = demuxEM_generate_gender_plot,
					sccloud_version = sccloud_version,
					zones = zones,
					num_cpu = num_cpu,
					memory = memory,
					disk_space = disk_space,
					preemptible = preemptible
			}
		}
	}

	if (generate_hashing_cite_seq_tasks.cite_seq_ids[0] != '') {
		scatter (cite_seq_id in generate_hashing_cite_seq_tasks.cite_seq_ids) {
			call tasks.run_scCloud_merge_rna_adt as run_scCloud_merge_rna_adt {
				input:
					input_raw_gene_bc_matrices_h5 = generate_hashing_cite_seq_tasks.id2rna[cite_seq_id],
					input_adt_csv = generate_hashing_cite_seq_tasks.id2adt[cite_seq_id],
					antibody_control_csv = antibody_control_csv,
					output_dir = sub(output_directory, "/+$", ""),
					output_name = cite_seq_id,
					sccloud_version = sccloud_version,
					zones = zones,
					memory = memory,
					disk_space = disk_space,
					preemptible = preemptible
			}
		}
	}
}
