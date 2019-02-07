workflow smartseq2_per_plate {
	# 3 columns (Cell, Read1, and Read2). gs URL
	File sample_sheet
	# Plate name
	String plate_name
	# Output directory, gs URL
	String output_directory

	# GRCh38, GRCm38 or a URL to a tar.gz file
	String reference

	File acronym_file = "gs://regev-lab/resources/SmartSeq2/index.tsv"
	# File acronym_file = "smartseq2_index.tsv"
	Map[String, String] acronym2gsurl = read_map(acronym_file)
	# If reference is a url
	Boolean is_url = sub(reference, "^.+\\.(tgz|gz)$", "URL") == "URL"

	File reference_file = (if is_url then reference else acronym2gsurl[reference])

	# smartseq2 version, default to "0.1.0"
	String? smartseq2_version = "0.1.0"
	# Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
	String? zones = "us-central1-b"
	# Number of cpus per job
	Int? num_cpu = 4
	# Memory in GB
	Int? memory = 10
	# disk space in GB
	Int? disk_space = 10
	# Number of preemptible tries 
	Int? preemptible = 2

	call parse_sample_sheet {
		input:
			sample_sheet = sample_sheet,
			smartseq2_version = smartseq2_version,
			zones = zones,
			preemptible = preemptible
	}

	scatter (i in range(length(parse_sample_sheet.cell_ids))) {
		call run_rsem {
			input:
				reference = reference_file,
				read1 = parse_sample_sheet.read1_list[i],
				read2 = parse_sample_sheet.read2_list[i],
				sample_name = parse_sample_sheet.cell_ids[i],
				smartseq2_version = smartseq2_version,
				zones = zones,
				num_cpu = num_cpu,
				memory = memory,
				disk_space = disk_space,
				preemptible = preemptible
		}
	}

	call generate_count_matrix {
		input:
			gene_results = run_rsem.rsem_gene,
			output_name = output_directory + "/" + plate_name,
			smartseq2_version = smartseq2_version,
			zones = zones,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible
	}
}

task parse_sample_sheet {
	File sample_sheet
	String smartseq2_version
	String zones
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import pandas as pd 
		from subprocess import check_call
		df = pd.read_csv('${sample_sheet}', header = 0, index_col = 0)
		with open('cell_ids.txt', 'w') as fo1, open('read1_list.txt', 'w') as fo2, open('read2_list.txt', 'w') as fo3:
			for cell_id, row in df.iterrows():
				fo1.write(cell_id + '\n')
				fo2.write(row['Read1'] + '\n')
				fo3.write(row['Read2'] + '\n')
		CODE
	}

	output {
		Array[String] cell_ids = read_lines('cell_ids.txt')
		Array[String] read1_list = read_lines('read1_list.txt')
		Array[String] read2_list = read_lines('read2_list.txt')
	}

	runtime {
		docker: "regevlab/smartseq2-${smartseq2_version}"
		zones: zones
		preemptible: "${preemptible}"
	}
}

task run_rsem {
	File reference
	File read1
	File read2
	String sample_name
	String smartseq2_version
	String zones	
	Int num_cpu
	Int memory
	Int disk_space
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &
		mkdir -p rsem_ref
		tar xf ${reference} -C rsem_ref --strip-components 1
		rsem-calculate-expression --bowtie2 --paired-end -p ${num_cpu} --append-names --time ${read1} ${read2} rsem_ref/rsem_ref ${sample_name}
	}

	output {
		File rsem_gene = "${sample_name}.genes.results"
		File rsem_isoform = "${sample_name}.isoforms.results"
		File rsem_time = "${sample_name}.time"
		File rsem_cnt = "${sample_name}.stat/${sample_name}.cnt"
		File rsem_model = "${sample_name}.stat/${sample_name}.model"
		File rsem_theta = "${sample_name}.stat/${sample_name}.theta"
	}

	runtime {
		docker: "regevlab/smartseq2-${smartseq2_version}"
		zones: zones
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: "${preemptible}"
	}
}

task generate_count_matrix {
	Array[File] gene_results
	String output_name
	String smartseq2_version
	String zones
	Int memory
	Int disk_space
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import os
		import numpy as np
		import pandas as pd
		gene_names = None
		barcodes = []
		cntmat = []
		for result_file in "${sep=',' gene_results}".split(','):
			barcodes.append(os.path.basename(result_file)[:-len('.genes.results')])
			df = pd.read_table(result_file, header = 0, index_col = 0)
			if gene_names is None:
				gene_names = np.array(['_'.join(x.split('_')[1:]) for x in df.index])
			tot_counts = df['expected_count'].sum()
			counts = df['TPM'].values.copy()
			counts = (counts / counts.sum() * tot_counts + 0.5).astype(int)
			cntmat.append(counts)
		df_idx = pd.Index(gene_names, name = 'GENE')
		df_out = pd.DataFrame(data = np.stack(cntmat, axis = 1), index = df_idx, columns = barcodes)
		df_out.to_csv('results.dge.txt.gz', sep = '\t', compression = 'gzip')
		CODE

		gsutil cp results.dge.txt.gz ${output_name}.dge.txt.gz
		# cp results.dge.txt.gz ${output_name}.dge.txt.gz
	}

	output {
		String output_count_matrix = "${output_name}.dge.txt.gz"
	}

	runtime {
		docker: "regevlab/smartseq2-${smartseq2_version}"
		zones: zones
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: "${preemptible}"
	}
}
