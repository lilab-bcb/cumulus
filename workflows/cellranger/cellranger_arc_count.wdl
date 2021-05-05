version 1.0

workflow cellranger_arc_count {
	input {
		# Sample ID
		String sample_id
		# A comma-separated list of input FASTQs directories (gs urls)
		String input_gex_fastqs_directories
		String input_atac_fastqs_directories

		String output_directory

		# Keywords or a URL to a tar.gz file
		String genome

		String cellranger_arc_version = "1.0.1"
		# Google cloud zones, default to "us-central1-b", which is consistent with Cromwell's genomics.default-zones attribute
		String zones = "us-central1-b"
		# Number of cpus per cellranger job
		Int num_cpu = 64
		# Memory string, e.g. 57.6G
		String memory = "57.6G"
		# Disk space in GB
		Int disk_space = 500
		# Number of preemptible tries
		Int preemptible = 2

		# Which docker registry to use: quay.io/cumulus (default) or cumulusprod
		String docker_registry = "quay.io/cumulus"
	}

	File acronym_file = "gs://regev-lab/resources/cellranger/index.tsv"
	# File acronym_file = "index.tsv"
	Map[String, String] acronym2gsurl = read_map(acronym_file)
	# If reference is a url
	Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"

	File genome_file = (if is_url then genome else acronym2gsurl[genome])

	call run_cellranger_arc_count {
		input:
			sample_id = sample_id,
			input_gex_fastqs_directories = input_gex_fastqs_directories,
			input_atac_fastqs_directories=input_atac_fastqs_directories,
			output_directory = sub(output_directory, "/+$", ""),
			genome_file = genome_file,
			cellranger_arc_version = cellranger_arc_version,
			zones = zones,
			num_cpu = num_cpu,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible,
			docker_registry = docker_registry
	}

	output {
		String output_count_directory = run_cellranger_arc_count.output_count_directory
		String output_metrics_summary = run_cellranger_arc_count.output_metrics_summary
		String output_web_summary = run_cellranger_arc_count.output_web_summary
	}
}

task run_cellranger_arc_count {
	input {
		String sample_id
		String input_gex_fastqs_directories
		String input_atac_fastqs_directories
		String output_directory
		File genome_file
		String cellranger_arc_version
		String zones
		Int num_cpu
		String memory
		Int disk_space
		Int preemptible
		String docker_registry
	}

	command <<<
		set -e

		monitor_script.sh &
		mkdir -p genome_dir
		tar xf ~{genome_file} -C genome_dir --strip-components 1

		python <<CODE
		import re
		import os
		from subprocess import check_call

		def download_fastqs(directories, offset):
			fastqs = []
			for i in range(len(directories)):
				directory = directories[i]
				directory = re.sub('/+$', '', directory) # remove trailing slashes
				dest = '~{sample_id}_' + str(i+offset)
				call_args = ['gsutil', '-m', 'cp', '-r', directory, '.']
				check_call(call_args)
				check_call(['mv', os.path.basename(directory), dest])
				fastqs.append(os.path.abspath(dest))
			return fastqs
		gex_fastqs = download_fastqs('~{input_gex_fastqs_directories}'.split(','), 0)
		atac_fastqs = download_fastqs('~{input_atac_fastqs_directories}'.split(','), len(gex_fastqs))
		with open('libraries.csv', 'wt') as out:
			out.write('fastqs,sample,library_type\n')
			for fq in gex_fastqs:
				out.write('{},{},{}\n'.format(fq, '~{sample_id}', 'Gene Expression'))
			for fq in atac_fastqs:
				out.write('{},{},{}\n'.format(fq, '~{sample_id}', 'Chromatin Accessibility'))

		call_args = ['cellranger-arc', 'count', '--id=results', '--libraries=$PWD/libraries.csv', '--reference=genome_dir', '--jobmode=local']
		check_call(call_args)
		CODE

		gsutil -q -m rsync -d -r results/outs "~{output_directory}"/~{sample_id}
		# cp -r results/outs "~{output_directory}"/~{sample_id}
	>>>

	output {
		String output_count_directory = "~{output_directory}/~{sample_id}"
		String output_metrics_summary = "~{output_directory}/~{sample_id}/summary.csv"
		String output_web_summary = "~{output_directory}/~{sample_id}/web_summary.html"
	}

	runtime {
		docker: "~{docker_registry}/cellranger-arc:~{cellranger_arc_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: "~{num_cpu}"
		preemptible: "~{preemptible}"
	}
}
