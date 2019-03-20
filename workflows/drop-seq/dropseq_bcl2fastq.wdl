workflow dropseq_bcl2fastq {
	String input_bcl_directory
	String output_directory

	# Whether to delete input bcl directory. If false, you should delete this folder yourself so as to not incur storage charges.
	Boolean? delete_input_bcl_directory = false
	String? zones = "us-central1-b"
	Int? num_cpu = 64
	String? memory = "128G"
	Int? disk_space = 1500
	Int? preemptible = 2
	Int minimum_trimmed_read_length = 10
	Int mask_short_adapter_reads = 10
	String? bcl2fastq_version = "2.20.0.422-2"

	call run_bcl2fastq {
		input:
			input_bcl_directory = sub(input_bcl_directory, "/+$", ""),
			output_directory = sub(output_directory, "/+$", ""),
			delete_input_bcl_directory = delete_input_bcl_directory,
			zones = zones,
			num_cpu = num_cpu,
			minimum_trimmed_read_length=minimum_trimmed_read_length,
			mask_short_adapter_reads=mask_short_adapter_reads,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible,
			bcl2fastq_version=bcl2fastq_version
	}

	output {
		String fastqs = run_bcl2fastq.fastqs
		File monitoringLog = run_bcl2fastq.monitoringLog
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
	Int minimum_trimmed_read_length
	Int mask_short_adapter_reads
	String run_id = basename(input_bcl_directory)
	String bcl2fastq_version


	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		gsutil -q -m cp -r ${input_bcl_directory} .

		cd ${run_id}

		bcl2fastq \
		--output-dir out \
		--no-lane-splitting \
		--minimum-trimmed-read-length ${minimum_trimmed_read_length} \
		--mask-short-adapter-reads ${mask_short_adapter_reads}

		cd out

		gsutil -q -m cp -r . ${output_directory}/${run_id}_fastqs/

		python <<CODE

		import os
		import re
		from subprocess import check_call
		gs_url = '${output_directory}/${run_id}_fastqs/'

		with open('../../sample_sheet.txt', 'w') as sample_sheet_writer:
            for project_dir in os.listdir('.'):
                if os.path.isdir(os.path.join(project_dir)) and project_dir != 'Reports' and project_dir != 'Stats':
                    for sample_id_dir in os.listdir(os.path.join(project_dir)):
                        fastq_files = os.listdir(os.path.join(project_dir, sample_id_dir))
                        if len(fastq_files) != 2:
                            raise ValueError(str(len(fastq_files)) + ' fastq files found')
                        fastq_files.sort()
                        fastq_path = gs_url + project_dir + '/' + sample_id_dir + '/'
                        sample_id = re.sub('_S[0-9]+_R1_001.fastq.gz', '',fastq_files[0])
                        sample_sheet_writer.write(sample_id)
                        sample_sheet_writer.write('\t' + fastq_path + fastq_files[0])
                        sample_sheet_writer.write('\t' + fastq_path + fastq_files[1])
                        sample_sheet_writer.write('\t' + str(os.path.getsize(fastq_files[0]) + os.path.getsize(fastq_files[1])))
                        sample_sheet_writer.write('\n')

		if '${delete_input_bcl_directory}' is 'true':
			call_args = ['gsutil', '-q', '-m', 'rm', '-r', '${input_bcl_directory}']
			check_call(call_args)
		CODE

		gsutil -q -m cp  ../../sample_sheet.txt "${output_directory}/${run_id}_fastqs.txt"
	}

	output {
		String fastqs = "${output_directory}/${run_id}_fastqs.txt"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/bcl2fastq-${bcl2fastq_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}
