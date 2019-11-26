workflow cellranger_atac_count {
	# Sample ID
	String sample_id
	# A comma-separated list of input FASTQs directories (gs urls)
	String input_fastqs_directories
	# cellRanger-atac output directory, gs url
	String output_directory

	# Keywords or a URL to a tar.gz file
	String genome


	File acronym_file = "gs://regev-lab/resources/cellranger/index.tsv"
	# File acronym_file = "index.tsv"
	Map[String, String] acronym2gsurl = read_map(acronym_file)
	# If reference is a url
	Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"

	File genome_file = (if is_url then genome else acronym2gsurl[genome])

	# Force pipeline to use this number of cells, bypassing the cell detection algorithm
	Int? force_cells
	# Chose the algorithm for dimensionality reduction prior to clustering and tsne: 'lsa' (default), 'plsa', or 'pca'.
	String? dim_reduce

	# 1.0.0 or 1.0.1
	String? cellranger_atac_version = "1.0.1"
	# Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
	String? zones = "us-central1-b"
	# Number of cpus per cellranger job
	Int? num_cpu = 64
	# Memory string, e.g. 57.6G
	String? memory = "57.6G"
	# Disk space in GB
	Int? disk_space = 500
	# Number of preemptible tries 
	Int? preemptible = 2
	String docker_registry

	call run_cellranger_atac_count {
		input:
			sample_id = sample_id,
			input_fastqs_directories = input_fastqs_directories,
			output_directory = sub(output_directory, "/+$", ""),
			genome_file = genome_file,
			force_cells=force_cells,
			dim_reduce = dim_reduce,
			cellranger_atac_version = cellranger_atac_version,
			zones = zones,
			num_cpu = num_cpu,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible,
			docker_registry = docker_registry
	}

	output {
		String output_count_directory = run_cellranger_atac_count.output_count_directory
		String output_metrics_summary = run_cellranger_atac_count.output_metrics_summary
		String output_web_summary = run_cellranger_atac_count.output_web_summary
		File monitoringLog = run_cellranger_atac_count.monitoringLog
	}
}

task run_cellranger_atac_count {
	String sample_id
	String input_fastqs_directories
	String output_directory
	File genome_file
	Int? force_cells
	String? dim_reduce
	String cellranger_atac_version
	String zones
	Int num_cpu
	String memory
	Int disk_space
	Int preemptible
	String docker_registry

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &
		mkdir -p genome_dir
		tar xf ${genome_file} -C genome_dir --strip-components 1

		python <<CODE
		import re
		from subprocess import check_call

		fastqs = []
		for i, directory in enumerate('${input_fastqs_directories}'.split(',')):
			directory = re.sub('/+$', '', directory) # remove trailing slashes
			call_args = ['gsutil', '-q', '-m', 'cp', '-r', directory + '/${sample_id}', '.']
			# call_args = ['cp', '-r', directory + '/${sample_id}', '.']
			print(' '.join(call_args))
			check_call(call_args)
			call_args = ['mv', '${sample_id}', '${sample_id}_' + str(i)]
			print(' '.join(call_args))
			check_call(call_args)
			fastqs.append('${sample_id}_' + str(i))

		call_args = ['cellranger-atac', 'count', '--id=results', '--reference=genome_dir', '--fastqs=' + ','.join(fastqs), '--sample=${sample_id}', '--jobmode=local']
		if '${force_cells}' is not '':
			call_args.append('--force-cells=${force_cells}')
		if '${dim_reduce}' is not '':
			call_args.append('--dim-reduce=${dim_reduce}')
		print(' '.join(call_args))
		check_call(call_args)
		CODE

		gsutil -q -m rsync -d -r results/outs ${output_directory}/${sample_id}
		# cp -r results/outs ${output_directory}/${sample_id}
	}

	output {
		String output_count_directory = "${output_directory}/${sample_id}"
		String output_metrics_summary = "${output_directory}/${sample_id}/summary.csv"
		String output_web_summary = "${output_directory}/${sample_id}/web_summary.html"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "${docker_registry}cellranger-atac:${cellranger_atac_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: "${preemptible}"
	}
}
