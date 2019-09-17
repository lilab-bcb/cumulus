workflow bcl2fastq {
	String input_bcl_directory
	String output_directory
	# Whether to delete input bcl directory. If false, you should delete this folder yourself so as to not incur storage charges.
	Boolean? delete_input_bcl_directory = false
	String? zones
	Int? num_cpu = 32
	String? memory = "120G"
	Int? disk_space = 1500
	Int? preemptible = 2
	Int? minimum_trimmed_read_length
	Int? mask_short_adapter_reads
	Int? barcode_mismatches
	Boolean? create_fastq_for_index_reads = false
	String? bcl2fastq_version = "2.20.0.422"
	String? docker_registry = "gcr.io/sccloud-prod"
	# If not specified, assume SampleSheet.csv is in input_bcl_directory
	File? sample_sheet

	call run_bcl2fastq {
		input:
			input_bcl_directory = sub(input_bcl_directory, "/+$", ""),
			output_directory = sub(output_directory, "/+$", ""),
			delete_input_bcl_directory = delete_input_bcl_directory,
			sample_sheet=sample_sheet,
			zones = zones,
			num_cpu = num_cpu,
			minimum_trimmed_read_length=minimum_trimmed_read_length,
			mask_short_adapter_reads=mask_short_adapter_reads,
			barcode_mismatches = barcode_mismatches,
			barcode_mismatches=barcode_mismatches,
			create_fastq_for_index_reads=create_fastq_for_index_reads,
			memory = memory,
			docker_registry = docker_registry,
			disk_space = disk_space,
			preemptible = preemptible,
			bcl2fastq_version=bcl2fastq_version
	}

	output {
		String fastqs = run_bcl2fastq.fastqs
	}
}

task run_bcl2fastq {
	String input_bcl_directory
	String output_directory
	Boolean delete_input_bcl_directory
	String zones
	Int num_cpu
	String memory
	Int disk_space
	Int preemptible
	Int? minimum_trimmed_read_length
	Int? mask_short_adapter_reads
	Int? barcode_mismatches
	Boolean create_fastq_for_index_reads
	String run_id = basename(input_bcl_directory)
	String bcl2fastq_version
	File? sample_sheet
	String docker_registry

	command {
		set -e
		export TMPDIR=/tmp

		gsutil -q -m cp -r ${input_bcl_directory} .

		cd ${run_id}

		bcl2fastq \
		--output-dir out \
		--no-lane-splitting \
		${"--sample-sheet " + sample_sheet} \
		${"--minimum-trimmed-read-length " + minimum_trimmed_read_length} \
		${"--mask-short-adapter-reads " + mask_short_adapter_reads} \
		${"--barcode-mismatches " + barcode_mismatches} \
		${true="--create-fastq-for-index-reads" false="" create_fastq_for_index_reads}

		cd out

		gsutil -q -m cp -r . ${output_directory}/${run_id}_fastqs/

		python <<CODE
		import os
		import re
		import subprocess

		gs_url = '${output_directory}/${run_id}_fastqs/'
		with open('../sample_sheet.txt', 'w') as sample_sheet_writer:
			sample_name_to_fastqs = dict()
			for root, dirs, files in os.walk("."):
				for name in files:
					if name.endswith('.fastq.gz') and not name.startswith('Undetermined'):
						sample_name = re.sub('_S[0-9]+_[R,I][1-9]+_001.fastq.gz$', '', name)
						fastqs = sample_name_to_fastqs.get(sample_name, None)
						if fastqs is None:
							fastqs = []
							sample_name_to_fastqs[sample_name] = fastqs
						fastq_path = os.path.normpath(os.path.join(root, name))
						fastqs.append(fastq_path)
			for sample_name in sample_name_to_fastqs:
				fastqs = sample_name_to_fastqs[sample_name]
				fastqs.sort()
				total_size = 0
				sample_sheet_writer.write(sample_name)
				for fastq in fastqs:
					total_size += os.path.getsize(fastq)
					sample_sheet_writer.write('\t' + gs_url + fastq)
				sample_sheet_writer.write('\t' + str(total_size) + '\n')

		if '${delete_input_bcl_directory}' is 'true':
			p = subprocess.run(['gsutil', '-q', '-m', 'rm', '-r', '${input_bcl_directory}'])
			if p.returncode != 0:
				print('Unable to delete BCL directory')

		CODE

		gsutil -q -m cp	 ../sample_sheet.txt "${output_directory}/${run_id}_fastqs.txt"
	}

	output {
		String fastqs = "${output_directory}/${run_id}_fastqs.txt"
	}

	runtime {
		docker: "${docker_registry}/bcl2fastq:${bcl2fastq_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}
