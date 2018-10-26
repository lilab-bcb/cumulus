import "https://api.firecloud.org/ga4gh/v1/tools/scCloud:tasks/versions/4/plain-WDL/descriptor" as tasks
# import "../scCloud/scCloud_tasks.wdl" as tasks

workflow scCloud_hashing_cite_seq {
	# An sample sheet contains RNA and ADT data correspondence
	File input_sample_sheet
	# Output directory, gs url
	String output_directory

	# Number of cpus
	Int? num_cpu = 8
	# Memory in GB
	Int? memory = 10
	# Disk space in GB
	Int? disk_space = 20
	# Number of preemptible tries 
	Int? preemptible = 2

	# Reference genome
	String? genome

	# demuxEM parameters
	# Only demultiplex cells/nuclei with at least <demuxEM_min_num_genes> expressed genes
	Int? demuxEM_min_num_genes
	# Any cell/nucleus with no less than <demuxEM_max_background_probability> background probability will be marked as unknown. [default: 0.8]
	Float? demuxEM_max_background_probability
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
					hash_type = generate_hashing_cite_seq_tasks.id2type[hashing_id],
					genome = genome,
					min_num_genes = demuxEM_min_num_genes,
					max_background_probability = demuxEM_max_background_probability,
					random_state = demuxEM_random_state,
					generate_diagnostic_plots = demuxEM_generate_diagnostic_plots,
					generate_gender_plot = demuxEM_generate_gender_plot,					
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
					memory = memory,
					disk_space = disk_space,
					preemptible = preemptible
			}
		}
	}
}
