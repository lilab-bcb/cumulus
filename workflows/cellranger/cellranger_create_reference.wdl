workflow cellranger_create_reference {
	String docker_registry = "cumulusprod/"
	String cellranger_version = '3.0.2'
	Int disk_space = 500
	Int preemptible = 2

	String input_file
	String output_dir
	String genome_list
	String fasta_list
	String genes_list
	String attributes

	Boolean is_sample_sheet = sub(input_file, "^.+\\.csv$", "CSV") == "CSV"


	if (is_sample_sheet)

	if ('${attributes}' != '') {
		call run_cellranger_filter as filter {
			input:
				docker_registry = docker_registry,
				cellranger_version = cellranger_version,
				disk_space = disk_space,
				preemptible = preemptible,
				input_gtf_file = intput_gtf_file,
				attributes = attributes
		}
	}




}

task generate_species {
	String docker_registry
	String cellranger_version
	String disk_space
	String preemptible
	String zones
	String memory

	String input_csv_file

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE

		import pandas as pd
		from subprocess import check_call

		df = pd.read_csv('${input_csv_file}', header = 0, dtype=str, index_col=False)
		for col in df.columns:
			df[col] = df[col].str.strip()



		CODE
	}


	runtime {
		docker: "${docker_registry}cellranger:${cellranger_version}"
		zones: zones
		preemptible: "${preemptible}"
	}
}

task run_cellranger_filter {
	String docker_registry
	String cellranger_version
	Int disk_space
	String zones
	String memory
	Int preemptible

	String input_gtf_file
	String output_name
	String attributes

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log

		python <<CODE
		from subprocess import check_call

		attrs = '${attributes}'.split(';')

		call_args = ['cellranger', 'mkgtf', '${input_gtf_file}', '${output_name}.filter.gtf']
		for attr in attrs:
			call_args.append('--attributes=' + attr)

		print(' '.join(call_args))
		check_all(call_args)
		CODE
	}

	output {
		File output_gtf_file = '${output_name}.filter.gtf'
	}

	runtime {
		docker: "${docker_registry}cellranger:${cellranger_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
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
	String memory
	Int preemptible

	String genome_list
	String fasta_list
	String genes_list
	String? ref_version

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log

		python <<CODE
		from subprocess import check_call

		call_args = ['cellranger', 'mkref']

		if '${ref_version}' is not '':
			call_args.append('--ref-version=${ref_version}')

		CODE
	}

	output {

	}

	runtime {
		docker: "${docker_registry}cellranger:${cellranger_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: "${preemptible}"
	}
}