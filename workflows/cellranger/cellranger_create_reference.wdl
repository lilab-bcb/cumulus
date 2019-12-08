workflow cellranger_create_reference {
	String? docker_registry = "cumulusprod/"
	String? cellranger_version = '3.1.0'
	Int? disk_space = 100
	Int? preemptible = 2
	String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	Int? num_cpu = 1
	Int? memory = 32

	File? input_sample_sheet
	File? input_gtf_file
	File? input_fasta
	String output_dir
	String? genome
	String? attributes
	String? ref_version

	if (input_sample_sheet != '') {
		call process_sample_sheet as pss {
			input:
				input_csv_file = input_sample_sheet
		}

		if (pss.genome_list[0] != '') {
			scatter (run_id in pss.genome_list) {
				call run_cellranger_filter as filt {
					input:
						docker_registry = docker_registry,
						cellranger_version = cellranger_version,
						disk_space = disk_space,
						zones = zones,
						memory = memory,
						preemptible = preemptible,
						input_gtf_file = pss.gtf_dict[run_id],
						attributes = pss.attr_dict[run_id],
						genome = run_id
				}
			}
		}

		call run_cellranger_create_reference_multi_species as crms {
			input:
				docker_registry = docker_registry,
				cellranger_version = cellranger_version,
				disk_space = disk_space,
				preemptible = preemptible,
				zones = zones,
				output_dir = output_dir,
				genomes = pss.genome_list,
				memory = memory,
				num_cpu = num_cpu,
				ref_version = ref_version
		}
	}

	if (input_sample_sheet == '') {
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
				fasta = input_fasta,
				genes = if do_filter then filter.output_gtf_file else input_gtf_file,
				memory = memory,
				num_cpu = num_cpu,
				ref_version = ref_version
		}
	}

}

task process_sample_sheet {
	File input_csv_file

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE

		import pandas as pd

		df = pd.read_csv('${input_csv_file}', header = 0, dtype = str, index_col=False)
		for c in df.columns:
			df[c] = df[c].str.strip()

		with open('genome_list'.txt, 'w') as fo1, open('fa_dict.txt', 'w') as fo2, open('gtf_dict.txt', 'w') as fo3, open('attr_dict.txt', 'w') as fo4:
			for _, row in df.iterrows():
				fo1.write(row['Genome'] + '\n')
				fo2.write(row['Genome'] + '\t' + row['Fasta'] + '\n')
				fo3.write(row['Genome'] + '\t' + row['Genes'] + '\n')
				fo4.write(row['Genome'] + '\t' + row['Attributes'] + '\n')

		CODE
	}

	output {
		Array[String] genome_list = read_lines('genome_list.txt')
		Map[String, String] fa_dict = read_map('fa_dict.txt')
		Map[String, String] gtf_dict = read_map('gtf_dict.txt')
		Map[String, String] attr_dict = read_map('attr_dict.txt')
	}
}

task run_cellranger_filter {
	String docker_registry
	String cellranger_version
	Int disk_space
	String zones
	Int memory
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

		do_filter = False
		attrs = '${attributes}'.split(';')

		call_args = ['cellranger', 'mkgtf', '${input_gtf_file}', '${genome}.filter.gtf']
		for attr in attrs:
			if attr is not '':
				call_args.append('--attribute=' + attr)
				if not do_filter: 
					do_filter = True

		if not do_filter:
			call_args = ['cp', '${input_gtf_file}', '${genome}.filter.gtf']
		
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

task run_cellranger_create_reference_multi_species {
	String docker_registry
	String cellranger_version
	Int disk_space
	Int preemptible
	String zones
	String output_dir
	Array[String] genomes
	Int memory
	Int num_cpu
	String ref_version

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call

		call_args = ['cellranger', 'mkref']
		
		genome_list = '${sep="," genomes}'.split(',')
		for g in genome_list:
			call_args.append('--genome=' + g)
			call_args.append('--fasta=' + g + '.fa')
			call_args.append('--genes=' + g + '.filter.gtf')

		call_args.append('--nthreads=${num_cpu}')
		call_args.append('--memgb=${memory}')

		if '${ref_version}' is not '':
			call_args.append('--ref-version=${ref_version}')

		print(' '.join(call_args))
		check_call(call_args)

		CODE

		tar -czf ${sep="_and_" genomes}.tar.gz ${sep="_and_" genomes}
		gsutil -q cp ${sep="_and_" genomes}.tar.gz ${output_dir}
		# mkdir -p ${output_dir}
		# cp ${sep="_and_" genomes}.tar.gz ${output_dir}

	}

	output {
		File reference = '${output_dir}/*.tar.gz'
		File monitoringLog = 'monitoring.log'
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
		File monitoringLog = "monitoring.log"
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