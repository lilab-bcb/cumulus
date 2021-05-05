version 1.0

workflow cellranger_multi {
	input {
		File csv
		String output_directory
		Int cpu = 30
		String memory = "45G"
		# Disk space in GB
		Int disk_space = 500
		Int preemptible = 2
		String docker = "cumulusprod/cellranger:6.0.0"
	}

	call run_cellranger_multi {
		input:
			csv = csv,
			acronym_file = "gs://regev-lab/resources/cellranger/index.tsv",
			output_directory = sub(output_directory, "/+$", ""),
			cpu = cpu,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible,
			docker = docker
	}

	output {
		String output_multi_directory = run_cellranger_multi.output_multi_directory
	}
}

task run_cellranger_multi {
	input {
		File csv
		String output_directory

		File acronym_file
		String memory
		Int disk_space
		Int cpu
		Int preemptible
		String docker
	}

	command <<<
		set -e

		monitor_script.sh &
		export PROJECT=$(gcloud config get-value project)

		python <<CODE
		import re
		import os
		from subprocess import check_call
		import pandas as pd

		project = os.getenv('PROJECT')
		csv = '~{csv}'
		output_directory = '~{output_directory}'
		acronym_file = '~{acronym_file}'
		acronym_df = pd.read_csv(acronym_file, header=None, index_col=0, sep='\t')

		def download_fastqs(sample_id, directory):
			directory = re.sub('/+$', '', directory) # remove trailing slashes
			offset = 1
			dest =	sample_id + '_' + str(offset)
			while os.path.exists(dest):
				offset += 1
				dest = sample_id + '_' + str(offset)
			call_args = ['gsutil', '-m', 'cp', '-r', directory, '.']
			check_call(call_args)
			check_call(['mv', os.path.basename(directory), dest])
			return os.path.abspath(dest)


		with open('new.csv', 'wt') as out:
			with open(csv, 'rt') as f:
				fastqs_index = -1
				fastq_id_index = -1
				section = None

				for line in f:
					line = line.strip()
					tokens = line.split(',')
					if len(tokens) > 0 and tokens[0].startswith('['):
						out.write(tokens[0] + '\n')
						section = tokens[0]
					else:
						if section == '[libraries]':
							if fastqs_index == -1:
								fastqs_index = tokens.index('fastqs')
								fastq_id_index = tokens.index('fastq_id')
							elif len(tokens) > 0 and tokens[0] != '':
								tokens[fastqs_index] = download_fastqs(tokens[fastq_id_index], tokens[fastqs_index])
						elif section == '[gene-expression]' and tokens[0] == 'reference':
							url_or_name = tokens[1]
							if url_or_name in acronym_df.index:
								url_or_name = acronym_df.loc[url_or_name][1]
							dest_tar = os.path.basename(url_or_name)
							check_call(['gsutil', '-u', project, '-m', 'cp', url_or_name, dest_tar])
							genome_dir = 'genome_dir'
							os.mkdir(genome_dir)
							check_call(['tar', 'xf', dest_tar, '-C', genome_dir, '--strip-components', '1'])
							tokens[1] = os.path.abspath('genome_dir')
						out.write(','.join(tokens) + '\n')

		call_args = ['cellranger', 'multi', '--id=results', '--csv=' + 'new.csv', '--jobmode=local']
		check_call(call_args)
		CODE

		gsutil -q -m rsync -d -r results/outs ~{output_directory}
	>>>

	output {
		String output_multi_directory = "~{output_directory}"
	}

	runtime {
		docker: docker
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: "~{cpu}"
		preemptible: "~{preemptible}"
	}
}
