workflow cellranger_create_reference {
	String? docker_registry = "cumulusprod/"
	String? cellranger_version = '3.1.0'
	Int? disk_space = 500
	Int? preemptible = 2
	String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	Int? num_cpu = 1
	Int? memory = 32

	File input_gtf_file
	String output_dir
	String genome
	File fasta
	String? attributes
	String? ref_version

	Boolean do_filter = if '${attributes}' != '' then true else false


	if (do_filter) {
		call run_cellranger_filter as filter {
			input:
				docker_registry = docker_registry,
				cellranger_version = cellranger_version,
				disk_space = disk_space,
				zones = zones,
				memory = memory,
				preemptible = preemptible,
				input_gtf_file = input_gtf_file,
				attributes = attributes,
				genome = genome
		}
	}


	call run_cellranger_create_reference as create_ref {
		input:
			docker_registry = docker_registry,
			cellranger_version = cellranger_version,
			disk_space = disk_space,
			preemptible = preemptible,
			zones = zones,
			output_dir = output_dir,
			genome = genome,
			fasta = fasta,
			genes = if do_filter then filter.output_gtf_file else input_gtf_file,
			memory = memory,
			num_cpu = num_cpu,
			ref_version = ref_version
	}

}

task run_cellranger_filter {
	String docker_registry
	String cellranger_version
	Int disk_space
	String zones
	String memory
	Int preemptible

	File input_gtf_file
	String genome
	String attributes

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call

		attrs = '${attributes}'.split(';')

		call_args = ['cellranger', 'mkgtf', '${input_gtf_file}', '${genome}.filter.gtf']
		for attr in attrs:
			call_args.append('--attribute=' + attr)

		print(' '.join(call_args))
		check_call(call_args)
		CODE
	}

	output {
		File output_gtf_file = "${genome}.filter.gtf"
	}

	runtime {
		docker: "${docker_registry}cellranger:${cellranger_version}"
		zones: zones
		memory: "${memory}G"
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: "${preemptible}"
	}
}

task run_cellranger_create_reference {
	String docker_registry
	String cellranger_version
	Int disk_space
	Int num_cpu
	String zones
	Int memory
	Int preemptible

	String output_dir
	String genome
	File fasta
	File genes
	String? ref_version

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call

		call_args = ['cellranger', 'mkref', '--genome=${genome}', '--fasta=${fasta}', '--genes=${genes}', '--nthreads=${num_cpu}', '--memgb=${memory}']

		if '${ref_version}' is not '':
			call_args.append('--ref-version=${ref_version}')

		print(' '.join(call_args))
		check_call(call_args)
		CODE

		tar -czf ${genome}.tar.gz ${genome}
		gsutil -q cp ${genome}.tar.gz ${output_dir}
		# mkdir -p ${output_dir}
		# cp ${genome}.tar.gz ${output_dir}
	}

	output {
		File reference = "${output_dir}/${genome}.tar.gz"
	}

	runtime {
		docker: "${docker_registry}cellranger:${cellranger_version}"
		zones: zones
		memory: "${memory}G"
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: "${preemptible}"
	}
}