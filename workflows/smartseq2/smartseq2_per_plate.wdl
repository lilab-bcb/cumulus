version 1.0

workflow smartseq2_per_plate {
	input {
		# 2-3 columns (Cell, Read1, and optionally Read2). gs URL
		File sample_sheet
		# Plate name
		String plate_name
		# Output directory, gs URL
		String output_directory

		# GRCh38, GRCm38 or a URL to a tar.gz file
		String reference
		
		# Align reads with 'aligner': hisat2-hca, star, bowtie2 (default hisat2-hca)
		String aligner = "hisat2-hca"

		# smartseq2 version
		String smartseq2_version = "1.1.0"
		# Google cloud zones, default to "us-central1-b", which is consistent with Cromwell's genomics.default-zones attribute
		String zones = "us-central1-b"
		# Number of cpus per job
		Int num_cpu = 4
		# Memory string
		String memory = "3.6G"
		# factor to multiply size of R1 and R2 by
	    Float disk_space_multiplier = 11
		# Disk space for count matrix generation task
	    Int generate_count_matrix_disk_space = 10
		# Number of preemptible tries 
		Int preemptible = 2
		# Which docker registry to use: cumulusprod (default) or quay.io/cumulus
	    String docker_registry = "cumulusprod"
	}


    File acronym_file = "gs://regev-lab/resources/smartseq2/index.tsv"
    # File acronym_file = "smartseq2_index.tsv"
    Map[String, String] acronym2gsurl = read_map(acronym_file)
    # If reference is a url
    Boolean is_url = sub(reference, "^.+\\.(tgz|gz)$", "URL") == "URL"

    String key = reference + "_" + aligner
    File reference_file = (if is_url then reference else acronym2gsurl[key])


	call parse_sample_sheet {
		input:
			sample_sheet = sample_sheet,
			smartseq2_version = smartseq2_version,
			zones = zones,
			preemptible = preemptible,
			docker_registry = docker_registry
	}

	# Paired-end data
	if (parse_sample_sheet.is_paired) {
		scatter (i in range(length(parse_sample_sheet.cell_ids))) {
			call run_rsem as run_rsem_pe {
				input:
					reference = reference_file,
					read1 = parse_sample_sheet.read1_list[i],
					read2 = parse_sample_sheet.read2_list[i],
					sample_name = parse_sample_sheet.cell_ids[i],
					aligner = aligner,
					smartseq2_version = smartseq2_version,
					zones = zones,
					num_cpu = num_cpu,
					memory = memory,
					disk_space_multiplier = disk_space_multiplier,
					preemptible = preemptible,
					docker_registry = docker_registry
			}
		}	
	}

	# Single-end data
	if (!parse_sample_sheet.is_paired) {
		scatter (i in range(length(parse_sample_sheet.cell_ids))) {
			call run_rsem as run_rsem_se {
				input:
					reference = reference_file,
					read1 = parse_sample_sheet.read1_list[i],
					sample_name = parse_sample_sheet.cell_ids[i],
					aligner = aligner,
					smartseq2_version = smartseq2_version,
					zones = zones,
					num_cpu = num_cpu,
					memory = memory,
					disk_space_multiplier = disk_space_multiplier,
					preemptible = preemptible,
					docker_registry = docker_registry
			}
		}		
	}

	call generate_count_matrix {
		input:
			gene_results = select_first([run_rsem_pe.rsem_gene, run_rsem_se.rsem_gene]),
			count_results = select_first([run_rsem_pe.rsem_cnt, run_rsem_se.rsem_cnt]),
			output_directory = output_directory,
			plate_name = plate_name,
			smartseq2_version = smartseq2_version,
			zones = zones,
			memory = memory,
			disk_space = generate_count_matrix_disk_space,
			preemptible = preemptible,
			docker_registry = docker_registry
	}

	output {
	    Array[File] rsem_gene = select_first([run_rsem_pe.rsem_gene, run_rsem_se.rsem_gene])
        Array[File] rsem_isoform = select_first([run_rsem_pe.rsem_isoform, run_rsem_se.rsem_isoform])
        Array[File] rsem_trans_bam = select_first([run_rsem_pe.rsem_trans_bam, run_rsem_se.rsem_trans_bam])
        Array[File] rsem_time = select_first([run_rsem_pe.rsem_time, run_rsem_se.rsem_time])
        Array[File] rsem_cnt = select_first([run_rsem_pe.rsem_cnt, run_rsem_se.rsem_cnt])
        Array[File] rsem_model = select_first([run_rsem_pe.rsem_model, run_rsem_se.rsem_model])
        Array[File] rsem_theta = select_first([run_rsem_pe.rsem_theta, run_rsem_se.rsem_theta])
		String output_count_matrix = generate_count_matrix.output_count_matrix
		String output_qc_report = generate_count_matrix.output_qc_report
	}
}


task parse_sample_sheet {
	input {
		File sample_sheet
		String smartseq2_version
		String zones
		Int preemptible
		String docker_registry	
	}

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import pandas as pd
		from subprocess import check_call
		df = pd.read_csv('~{sample_sheet}', header = 0, index_col = 0)
		is_paired = 'Read2' in df.columns
		with open('is_paired.txt', 'w') as fo:
			fo.write('true' if is_paired else 'false')
		with open('cell_ids.txt', 'w') as fo1, open('read1_list.txt', 'w') as fo2, open('read2_list.txt', 'w') as fo3:
			for cell_id, row in df.iterrows():
				fo1.write(cell_id + '\n')
				fo2.write(row['Read1'] + '\n')
				if is_paired:
					fo3.write(row['Read2'])
				fo3.write('\n')
		CODE
	}

	output {
		Array[String] cell_ids = read_lines('cell_ids.txt')
		Array[String] read1_list = read_lines('read1_list.txt')
		Array[String] read2_list = read_lines('read2_list.txt')
		Boolean is_paired = read_boolean('is_paired.txt')
	}

	runtime {
		docker: "~{docker_registry}/smartseq2:~{smartseq2_version}"
		zones: zones
		preemptible: preemptible
	}
}

task run_rsem {
	input {
		File reference
		File read1
		File? read2
		String sample_name
		String aligner
		String smartseq2_version
		String zones
		Int num_cpu
		String memory
		Float disk_space_multiplier
		Int preemptible
		String docker_registry	
	}

	Boolean is_star = aligner == "star"

	command {
		set -e
		export TMPDIR=/tmp

		mkdir -p rsem_ref
		tar xf ${reference} -C rsem_ref --strip-components 1
		REFERENCE_NAME="$(basename `ls rsem_ref/*.grp` .grp)"
		echo $REFERENCE_NAME
		rsem-calculate-expression --~{aligner} ~{true="--star-gzipped-read-file" false="" is_star} ~{true="--paired-end" false="" defined(read2)} -p ~{num_cpu} --append-names --time ~{read1} ~{default="" read2} rsem_ref/$REFERENCE_NAME ~{sample_name}
	}

	output {
		File rsem_gene = "~{sample_name}.genes.results"
		File rsem_isoform = "~{sample_name}.isoforms.results"
		File rsem_trans_bam = "~{sample_name}.transcript.bam"
		File rsem_time = "~{sample_name}.time"
		File rsem_cnt = "~{sample_name}.stat/~{sample_name}.cnt"
		File rsem_model = "~{sample_name}.stat/~{sample_name}.model"
		File rsem_theta = "~{sample_name}.stat/~{sample_name}.theta"
	}

	runtime {
		docker: "~{docker_registry}/smartseq2:~{smartseq2_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk " + ceil(size(reference, "GB")*5 + (disk_space_multiplier * (size(read1, "GB") + size(read2, "GB"))) + 1)+ " HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task generate_count_matrix {
	input {
		Array[File] gene_results
		Array[File] count_results
		String output_directory
		String plate_name
		String smartseq2_version
		String zones
		String memory
		Int disk_space
		Int preemptible
		String docker_registry	
	}

	String output_name = output_directory + "/" + plate_name

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
		for result_file in "~{sep=',' gene_results}".split(','):
			barcodes.append(os.path.basename(result_file)[:-len('.genes.results')])
			df = pd.read_table(result_file, header = 0, index_col = 0)
			if gene_names is None:
				gene_names = np.array(['_'.join(x.split('_')[1:]) for x in df.index])
			tot_counts = df['expected_count'].sum()
			counts = df['TPM'].values / 10.0 # convert TPMs into TP100Ks
			denom = counts.sum()
			if denom > 0:
				counts = (counts / denom * tot_counts + 0.5).astype(int)
			cntmat.append(counts)
		df_idx = pd.Index(gene_names, name = 'GENE')
		df_out = pd.DataFrame(data = np.stack(cntmat, axis = 1), index = df_idx, columns = barcodes)
		df_out.to_csv('results.dge.txt.gz', sep = '\t', compression = 'gzip')

		arr = []
		barcodes = []
		for result_file in "~{sep=',' count_results}".split(','):
			barcodes.append(os.path.basename(result_file)[:-len('.cnt')])
			with open(result_file) as fin:
				Ns = [int(x) for x in next(fin).strip().split(' ')]
				align_values = [int(x) for x in next(fin).strip().split(' ')]
				res = [str(Ns[3]), str(round(Ns[1] * 100.0 / Ns[3], 2)) + "%", str(round(align_values[0] * 100.0 / Ns[3], 2)) + "%"]
				arr.append(res)
		df = pd.DataFrame(data = np.array(arr), index = barcodes, columns = ["Total reads", "Alignment rate", "Unique rate"])
		df.index.name = "Cell"
		df.to_csv('results.qc.stat.tsv', sep = '\t')
		CODE

		gsutil cp results.dge.txt.gz ~{output_name}.dge.txt.gz
		gsutil cp results.qc.stat.tsv ~{output_name}.qc.stat.tsv
		# mkdir -p ~{output_directory}
		# cp results.dge.txt.gz ~{output_name}.dge.txt.gz
		# cp results.qc.stat.tsv ~{output_name}.qc.stat.tsv
	}

	output {
		String output_count_matrix = "~{output_name}.dge.txt.gz"
		String output_qc_report = "~{output_name}.qc.stat.tsv"
	}

	runtime {
		docker: "~{docker_registry}/smartseq2:~{smartseq2_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}
