workflow scCloud_adt {
	# Sample ID
	String sample_id
	# A comma-separated list of input FASTQs directories (gs urls)
	String input_fastqs_directories
	# Output directory, gs url
	String output_directory

	# 10x genomics chemistry 
	String chemistry

	# data type, either adt or crispr
	String data_type

	# cell barcodes white list, from 10x genomics, can be either v2 or v3 chemistry
	File cell_barcode_file = (if chemistry == "SC3Pv3" then "gs://regev-lab/resources/cellranger/3M-february-2018.txt.gz" else "gs://regev-lab/resources/cellranger/737K-august-2016.txt.gz")
	# File cell_barcode_file = (if chemistry == "SC3Pv3" then "3M-february-2018.txt.gz" else "737K-august-2016.txt.gz")

	# feature barcodes in csv format
	File feature_barcode_file

	# maximum hamming distance in feature barcodes
	Int? max_mismatch = 3

	# Memory in GB
	Int? memory = 32
	# Disk space in GB
	Int? disk_space = 100
	# Number of preemptible tries 
	Int? preemptible = 2

	call run_generate_count_matrix_ADTs {
		input:
			sample_id = sample_id,
			input_fastqs_directories = input_fastqs_directories,
			output_directory = sub(output_directory, "/+$", ""),
			chemistry = chemistry,
			data_type = data_type,
			cell_barcodes = cell_barcode_file,
			feature_barcodes = feature_barcode_file,
			max_mismatch = max_mismatch,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible
	}

	output {
		String output_count_directory = run_generate_count_matrix_ADTs.output_count_directory
		File monitoringLog = run_generate_count_matrix_ADTs.monitoringLog
	}
}

task run_generate_count_matrix_ADTs {
	String sample_id
	String input_fastqs_directories
	String output_directory
	String chemistry
	String data_type
	File cell_barcodes
	File feature_barcodes
	Int max_mismatch
	Int memory
	Int disk_space
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

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
	
		call_args = ['generate_count_matrix_ADTs', '${cell_barcodes}', '${feature_barcodes}', ','.join(fastqs), '${sample_id}', '--max-mismatch-feature', '${max_mismatch}']
		if '${data_type}' is 'crispr':
			call_args.extend(['--feature', 'crispr', '--min-reads-per-umi', '11', '--min-ratio-per-umi', '0.25'])
		else:
			call_args.extend(['--feature', 'antibody'])
		if '${chemistry}' is 'SC3Pv3':
			call_args.extend(['--max-mismatch-cell', '0', '--umi-length', '12'])
		else:
			call_args.extend(['--max-mismatch-cell', '1', '--umi-length', '10'])
		print(' '.join(call_args))
		check_call(call_args)
		CODE

		gsutil -q -m cp ${sample_id}.*csv* ${output_directory}/${sample_id}/
		# mkdir -p ${output_directory}/${sample_id}
		# cp -f ${sample_id}.*csv* ${output_directory}/${sample_id}/
	}

	output {
		String output_count_directory = "${output_directory}/${sample_id}"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/sccloud-0.6.0"
		zones: "us-central1-c us-central1-b us-east1-b us-east1-c us-east1-d"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: "${preemptible}"
	}
}
