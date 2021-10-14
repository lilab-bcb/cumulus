version 1.0

workflow cumulus_adt {
	input {
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

		# feature barcodes in csv format
		File feature_barcode_file

		# Index TSV file
		File acronym_file

		# scaffold sequence for Perturb-seq, default is "", which for Perturb-seq means barcode starts at position 0 of read 2
		String scaffold_sequence = ""

		# maximum hamming distance in feature barcodes
		Int max_mismatch = 3
		# minimum read count ratio (non-inclusive) to justify a feature given a cell barcode and feature combination, only used for crispr
		Float min_read_ratio = 0.1

		# cumulus_feature_barcoding version
		String cumulus_feature_barcoding_version
		# Which docker registry to use: quay.io/cumulus (default), or cumulusprod
		String docker_registry = "quay.io/cumulus"

		# Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
		String zones = "us-central1-b"
		# Memory string, e.g. 32G
		String memory = "32G"
		# Disk space in GB
		Int disk_space = 100

		# Number of preemptible tries
		Int preemptible = 2
		# Number of maximum retries when running on AWS
		Int awsMaxRetries = 5
		# Backend
		String backend = "gcp"
	}

	call run_generate_count_matrix_ADTs {
		input:
			sample_id = sample_id,
			input_fastqs_directories = input_fastqs_directories,
			output_directory = output_directory,
			chemistry = chemistry,
			data_type = data_type,
			acronym_file = acronym_file,
			feature_barcodes = feature_barcode_file,
			scaffold_sequence = scaffold_sequence,
			max_mismatch = max_mismatch,
			min_read_ratio = min_read_ratio,
			cumulus_feature_barcoding_version = cumulus_feature_barcoding_version,
			docker_registry = docker_registry,
			zones = zones,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible,
			awsMaxRetries = awsMaxRetries,
			backend = backend
	}

	output {
		String output_count_directory = run_generate_count_matrix_ADTs.output_count_directory
		File monitoringLog = run_generate_count_matrix_ADTs.monitoringLog
	}
}

task run_generate_count_matrix_ADTs {
	input {
		String sample_id
		String input_fastqs_directories
		String output_directory
		String chemistry
		String data_type
		File feature_barcodes
		File acronym_file
		String scaffold_sequence
		Int max_mismatch
		Float min_read_ratio
		String cumulus_feature_barcoding_version
		String docker_registry
		String zones
		String memory
		Int disk_space
		Int preemptible
		Int awsMaxRetries
		String backend
	}

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		import re
		from subprocess import check_call
		import pandas as pd

		fastqs = []

		ref_data = pd.read_csv('~{acronym_file}', sep = '\t', header = None, names = ['Keyword', 'Path'])
		idx = ref_data[ref_data['Keyword'] == '~{chemistry}'].index.values
		cell_barcodes = ref_data.at[idx[0],'Path']

		for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
			directory = re.sub('/+$', '', directory) # remove trailing slashes
			call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', '-r', directory + '/~{sample_id}/', './~{sample_id}/']
			print(' '.join(call_args))
			check_call(call_args)
			call_args = ['mv', '~{sample_id}', '~{sample_id}_' + str(i)]
			print(' '.join(call_args))
			check_call(call_args)
			fastqs.append('~{sample_id}_' + str(i))

                call_args = ['strato', 'cp', '--backend', '-m', '~{backend}', cell_barcodes, '.']
                print(' '.join(call_args))
                check_call(call_args)
                
                cell_barcodes_file = os.path.basename(cell_barcodes)
		call_args = ['generate_count_matrix_ADTs', cell_barcodes_file, '~{feature_barcodes}', ','.join(fastqs), '~{sample_id}', '--max-mismatch-feature', '~{max_mismatch}']
		if '~{data_type}' == 'crispr':
			call_args.extend(['--feature', 'crispr', '--scaffold-sequence', '~{scaffold_sequence}'])
			if '~{chemistry}' != 'SC3Pv3':
				call_args.append('--no-match-tso')
		else:
			call_args.extend(['--feature', 'antibody'])
			if '~{data_type}' == 'cmo':
				call_args.append('--convert-cell-barcode')
		if '~{chemistry}' == 'SC3Pv3':
			call_args.extend(['--max-mismatch-cell', '0', '--umi-length', '12'])
		else:
			call_args.extend(['--max-mismatch-cell', '1', '--umi-length', '10'])
		print(' '.join(call_args))
		check_call(call_args)
		CODE

		if [ -f "~{sample_id}".stat.csv.gz ]
		then
			filter_chimeric_reads ~{data_type} ~{feature_barcodes} "~{sample_id}.stat.csv.gz" ~{min_read_ratio} ~{sample_id}
		fi

		#gsutil -m cp "~{sample_id}".*csv* "~{output_directory}/~{sample_id}/"
		strato cp --backend ~{backend} -m "~{sample_id}".*csv* "~{output_directory}/~{sample_id}/"
		# mkdir -p "~{output_directory}/~{sample_id}"
		# cp -f "~{sample_id}".*csv* "~{output_directory}/~{sample_id}/"

		if [ -f "~{sample_id}".umi_count.pdf ]
		then
			# gsutil cp "~{sample_id}".umi_count.pdf "~{output_directory}/~{sample_id}/"
			strato cp --backend ~{backend} "~{sample_id}".umi_count.pdf "~{output_directory}/~{sample_id}/"
			# cp -f "~{sample_id}".umi_count.pdf "~{output_directory}/~{sample_id}/"
		fi
	}

	output {
		String output_count_directory = "~{output_directory}/~{sample_id}"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "~{docker_registry}/cumulus_feature_barcoding:~{cumulus_feature_barcoding_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: 1
		preemptible: preemptible
		maxRetries: if backend == "aws" then awsMaxRetries else 0
	}
}
