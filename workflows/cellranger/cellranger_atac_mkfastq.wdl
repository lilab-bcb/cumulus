workflow cellranger_atac_mkfastq {
	# Input BCL directory, gs url
	String input_bcl_directory
	# 3 column CSV file (Lane, Sample, Index)
	File input_csv_file
	# CellRanger output directory, gs url
	String output_directory

	# Whether to delete input bcl directory. If false, you should delete this folder yourself so as to not incur storage charges.
	Boolean? delete_input_bcl_directory = true
	# 1.0.0 or 1.0.1
	String? cellranger_atac_version = "1.0.1"
	# Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
	String? zones = "us-central1-b"
	# Number of cpus per cellranger-atac job
	Int? num_cpu = 64
	# Memory string, e.g. 128G
	String? memory = "128G"
	# Disk space in GB
	Int? disk_space = 1500
	# Number of preemptible tries 
	Int? preemptible = 2

	call run_cellranger_atac_mkfastq {
		input:
			input_bcl_directory = sub(input_bcl_directory, "/+$", ""),
			input_csv_file = input_csv_file,
			output_directory = sub(output_directory, "/+$", ""),
			delete_input_bcl_directory = delete_input_bcl_directory,
			cellranger_atac_version = cellranger_atac_version,
			zones = zones,
			num_cpu = num_cpu,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible
	}

	output {
		String output_fastqs_directory = run_cellranger_atac_mkfastq.output_fastqs_directory
		String output_fastqs_flowcell_directory = run_cellranger_atac_mkfastq.output_fastqs_flowcell_directory
		File monitoringLog = run_cellranger_atac_mkfastq.monitoringLog
	}
}

task run_cellranger_atac_mkfastq {
	String input_bcl_directory
	File input_csv_file
	String output_directory
	Boolean delete_input_bcl_directory
	String cellranger_atac_version
	String zones
	Int num_cpu
	String memory
	Int disk_space
	Int preemptible

	String run_id = basename(input_bcl_directory)

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &
		gsutil -q -m cp -r ${input_bcl_directory} .
		# cp -r ${input_bcl_directory} .
		cellranger-atac mkfastq --id=results --run=${run_id} --csv=${input_csv_file} --jobmode=local --qc

		python <<CODE
		import os
		import glob
		import pandas as pd
		from subprocess import check_call
		with open("output_fastqs_flowcell_directory.txt", "w") as fout:
			flowcell = [name for name in os.listdir('results/outs/fastq_path') if name != 'Reports' and name != 'Stats' and os.path.isdir('results/outs/fastq_path/' + name)][0]
			fout.write('${output_directory}/${run_id}_fastqs/fastq_path/' + flowcell + '\n')
		CODE

		gsutil -q -m rsync -d -r results/outs ${output_directory}/${run_id}_fastqs
		# cp -r results/outs ${output_directory}/${run_id}_fastqs

		python <<CODE
		from subprocess import check_call, check_output, CalledProcessError
		if '${delete_input_bcl_directory}' is 'true':
			try:
				call_args = ['gsutil', '-q', 'stat', '${output_directory}/${run_id}_fastqs/input_samplesheet.csv']
				print(' '.join(call_args))
				check_output(call_args)
				call_args = ['gsutil', '-q', '-m', 'rm', '-r', '${input_bcl_directory}']
				print(' '.join(call_args))
				check_call(call_args)
				print('${input_bcl_directory} is deleted!')
			except CalledProcessError:
				print("Failed to move outputs to Google bucket.")
		CODE
	}

	output {
		String output_fastqs_directory = "${output_directory}/${run_id}_fastqs"
		String output_fastqs_flowcell_directory = select_first(read_lines("output_fastqs_flowcell_directory.txt"))
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/cellranger-atac-${cellranger_atac_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible		
	}
}
