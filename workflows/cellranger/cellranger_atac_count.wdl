workflow cellranger_atac_count {
	File acronym_file = "gs://regev-lab/resources/cellranger-atac/index.tsv"

   # mm10v1.0.1 or a URL to a tar.gz file
	String genome
	Map[String, String] acronym2gsurl = read_map(acronym_file)
		# If reference is a url
	Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"

	File genome_file = (if is_url then genome else acronym2gsurl[genome])
	String output_directory
	String input_fastqs_directories
	String sample_id

	# Force pipeline to use this number of cells, bypassing the cell detection algorithm
	Int? force_cells
	String? cellranger_version = "1.0.1"

	# Number of cpus per cellranger job
	Int? num_cpu = 32
	# Memory in GB
	Int? memory = 64
	# Disk space in GB
	Int? disk_space = 500
	# Number of preemptible tries
	Int? preemptible = 2

	call run_cellranger_count {
		input:
			sample_id=sample_id,
			input_fastqs_directories=input_fastqs_directories,
			output_directory=output_directory,
			genome_file=genome_file,
			force_cells=force_cells,
			preemptible=preemptible,
			cpu=num_cpu,
			memory=memory,
			disk_space=disk_space,
			cellranger_version=cellranger_version
	}

}

task run_cellranger_count {
	File genome_file
	String sample_id
	String input_fastqs_directories
	String output_directory

	Int? force_cells
	Int cpu
	Int memory
	Int disk_space
	Int preemptible
	String cellranger_version

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
			check_call(call_args)
			call_args = ['mv', '${sample_id}', '${sample_id}_' + str(i)]
			check_call(call_args)
			fastqs.append('${sample_id}_' + str(i))

		call_args = ['cellranger-atac', 'count', '--id=results', '--reference=genome_dir', '--fastqs=' + ','.join(fastqs), '--sample=${sample_id}',  '--jobmode=local']
		if '${force_cells}' is not '':
			call_args.append('--force-cells=${force_cells}')

		check_call(call_args)
		CODE

		gsutil -q -m rsync -d -r results/outs ${output_directory}/${sample_id}
		# cp -r results/outs ${output_directory}/${sample_id}
	}

	output {
		String output_count_directory = "${output_directory}/${sample_id}"
		String output_metrics_summary = "${output_directory}/${sample_id}/metrics_summary.csv"
		String output_web_summary = "${output_directory}/${sample_id}/web_summary.html"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/cellranger-atac-${cellranger_version}"
		memory: "${memory} GB"
		disks: "local-disk ${disk_space} HDD"
		cpu: "${cpu}"
		preemptible: "${preemptible}"
	}
}
