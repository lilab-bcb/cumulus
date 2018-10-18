workflow cellranger_vdj {
	# Sample ID
	String sample_id
	# A comma-separated list of input FASTQs directories (gs urls)
	String input_fastqs_directories
	# CellRanger output directory, gs url
	String output_directory

	# GRCh38_vdj, GRCm38_vdj or a URL to a tar.gz file
	String genome


	File acronym_file = "gs://regev-lab/resources/cellranger/index.tsv"
	# File acronym_file = "index.tsv"
	Map[String, String] acronym2gsurl = read_map(acronym_file)
	# If reference is a url
	Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"

	File genome_file = (if is_url then genome else acronym2gsurl[genome])

	# Force pipeline to use this number of cells, bypassing the cell detection algorithm.
	Int? force_cells
	# Do not align reads to reference V(D)J sequences before de novo assembly. Default: false
	Boolean? denovo = false
	# Force the web summary HTML and metrics summary CSV to only report on a particular chain type. The accepted values are: auto for autodetection based on TR vs IG representation, TR for T cell receptors, IG for B cell receptors, all for all chain types.
	String? chain 

	# 2.2.0
	String? cellranger_version = "2.2.0"

	# Number of cpus per cellranger job
	Int? num_cpu = 64
	# Memory in GB
	Int? memory = 128
	# Disk space in GB
	Int? disk_space = 500
	# Number of preemptible tries 
	Int? preemptible = 2

	call run_cellranger_vdj {
		input:
			sample_id = sample_id,
			input_fastqs_directories = input_fastqs_directories,
			output_directory = sub(output_directory, "/+$", ""),
			genome_file = genome_file,
			force_cells = force_cells,
			denovo = denovo,
			chain = chain,
			cellranger_version = cellranger_version,
			num_cpu = num_cpu,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible
	}

	output {
		String output_vdj_directory = run_cellranger_vdj.output_vdj_directory
		String output_metrics_summary = run_cellranger_vdj.output_metrics_summary
		String output_web_summary = run_cellranger_vdj.output_web_summary
		File monitoringLog = run_cellranger_vdj.monitoringLog
	}
}

task run_cellranger_vdj {
	String sample_id
	String input_fastqs_directories
	String output_directory
	File genome_file
	Int? force_cells
	Boolean denovo
	String? chain
	String cellranger_version
	Int num_cpu
	Int memory
	Int disk_space
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &
		mkdir -p ref_dir
		tar xf ${genome_file} -C ref_dir --strip-components 1

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

		call_args = ['cellranger', 'vdj', '--id=results', '--reference=ref_dir', '--fastqs=' + ','.join(fastqs), '--sample=${sample_id}', '--jobmode=local']
		if '${force_cells}' is not '':
			call_args.append('--force-cells=${force_cells}')
		if '${denovo}' is not 'false':
			call_args.append('--denovo')
		if '${chain}' in ['auto', 'TR', 'IG', 'all']:
			call_args.append('--chain=${chain}')
		elif '${chain}' is not '':
			print('Unrecognized --chain value: ${chain}!')
		print(' '.join(call_args))
		check_call(call_args)
		CODE

		gsutil -q -m rsync -d -r results/outs ${output_directory}/${sample_id}
		# cp -r results/outs ${output_directory}/${sample_id}
	}

	output {
		String output_vdj_directory = "${output_directory}/${sample_id}"
		String output_metrics_summary = "${output_directory}/${sample_id}/metrics_summary.csv"
		String output_web_summary = "${output_directory}/${sample_id}/web_summary.html"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: "${preemptible}"
	}
}
