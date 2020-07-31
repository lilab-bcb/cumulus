version 1.0

workflow bulk_rna_seq {
	input {
		File sample_sheet
		String reference
		Boolean output_genome_bam = false
		Int num_cpu = 4
		String memory = "32G"
		# factor to multiply size of R1 and R2 by
		Float disk_space_multiplier = 4
		Float extra_disk_space = 2
		Int aggregate_extra_disk_space = 1
		Int preemptible = 2

		String rsem_docker = "cumulusprod/bulk_rna_seq:1.0.0"
		String util_docker = "cumulusprod/bulk_rna_seq_util:1.0.0"
		Map[String, String] qc_vars = {"mitochrondrial":"MT-", "ribosome":"RPL,RPS"}
	}


	File acronym_file = "gs://regev-lab/resources/bulk_rna_seq/index.json"

    Object acronym2gsurl = read_json(acronym_file)
	# If reference is a url
	Boolean is_url = sub(reference, "^.+\\.(tgz|gz)$", "URL") == "URL"
	File reference_file = (if is_url then reference else acronym2gsurl[reference]["star_genome"])

	call parse_sample_sheet {
		input:
			sample_sheet = sample_sheet,
			docker = util_docker,
			preemptible = preemptible
	}

	 scatter (name in parse_sample_sheet.names) {
	    call run_rsem  {
            input:
                reference = reference_file,
                read1 = parse_sample_sheet.r1[name],
                read2 = parse_sample_sheet.r2[name],
                extra_disk_space = extra_disk_space,
                sample_name = name,
                aligner = "star",
                output_genome_bam = output_genome_bam,
                num_cpu = num_cpu,
                memory = memory,
                disk_space_multiplier = disk_space_multiplier,
                preemptible = preemptible,
                docker = rsem_docker
        }
	}


	call aggregate {
		input:
			gene_results = run_rsem.rsem_gene,
			count_results = run_rsem.rsem_tpm,
			qc_vars = qc_vars,
			docker = util_docker,
			memory = memory,
			extra_disk_space = aggregate_extra_disk_space,
			preemptible = preemptible
	}


	output {
		Array[File] rsem_gene = run_rsem.rsem_gene
		Array[Array[File]] rsem_genome_bam  = run_rsem.rsem_genome_bam
		Array[File] rsem_isoform = run_rsem.rsem_isoform
		Array[File] rsem_trans_bam = run_rsem.rsem_trans_bam
		Array[File] rsem_time = run_rsem.rsem_time
		Array[File] aligner_log = run_rsem.aligner_log
		Array[File] rsem_tpm = run_rsem.rsem_tpm
		Array[File] rsem_model = run_rsem.rsem_model
		Array[File] rsem_theta = run_rsem.rsem_theta
		File rsem_agg_expected_counts = aggregate.rsem_agg_expected_counts
		File stats = aggregate.stats
	}
}

task parse_sample_sheet {
	input {
		File sample_sheet
		Int preemptible
		String docker
	}

	command {
		set -e
		sample_sheet.py ~{sample_sheet}
	}

	output {
		Array[String] names = read_lines('names.txt')
		Map[String, Array[String]] r1 = read_json('r1.json')
		Map[String, Array[String]] r2 = read_json('r2.json')
	}

	runtime {
		docker: docker
		preemptible: preemptible
	}
}

task run_rsem {
	input {
		File reference
		Array[File] read1
		Array[File] read2
		Boolean output_genome_bam
		String sample_name
		String aligner
		Float extra_disk_space
		Int num_cpu
		String memory
		Float disk_space_multiplier
		Int preemptible
		String docker
	}

	Boolean is_star = aligner == "star"
	Boolean paired_end = length(read2) > 0
	Boolean is_gzipped = sub(read1[0], "^.+\\.(gz)$", "GZ") == "GZ"
	Boolean star_gzipped_read_file = is_star && is_gzipped

	command {
		set -e

		monitor_script.sh > monitoring.log &

		mkdir -p rsem_ref
		tar xf ~{reference} -C rsem_ref --strip-components 1
		REFERENCE_NAME="$(basename `ls rsem_ref/*.grp` .grp)"
		rsem-calculate-expression --~{aligner} ~{true="--output-genome-bam" false="" output_genome_bam} ~{true="--star-gzipped-read-file" false="" star_gzipped_read_file} ~{true="--paired-end" false="" paired_end} -p ~{num_cpu} --append-names --time ~{sep=',' read1} ~{sep=',' read2} rsem_ref/$REFERENCE_NAME ~{sample_name}
	}

	output {
		File rsem_gene = "~{sample_name}.genes.results"
		File rsem_isoform = "~{sample_name}.isoforms.results"
		File rsem_trans_bam = "~{sample_name}.transcript.bam"
		File rsem_time = "~{sample_name}.time"
		File aligner_log = "~{sample_name}.log"
		File rsem_tpm = "~{sample_name}.stat/~{sample_name}.cnt"
		File rsem_model = "~{sample_name}.stat/~{sample_name}.model"
		File rsem_theta = "~{sample_name}.stat/~{sample_name}.theta"
		Array[File] rsem_genome_bam = glob("~{sample_name}.genome.bam")
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: docker
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk " + ceil(size(reference, "GB")*5 + (disk_space_multiplier * (size(read1, "GB") + size(read2, "GB"))) + extra_disk_space)+ " HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}


task aggregate {
	input {
		Array[File] gene_results
		Array[File] count_results
		Map[String, String] qc_vars
		String memory
		Int extra_disk_space
		Int preemptible
		String docker
	}


	command {
		set -e

		rsem_agg.py --gene ~{sep=' ' gene_results} --count ~{sep=' ' count_results} --qc_vars ~{write_json(qc_vars)}
	}

	output {
		File rsem_agg_expected_counts = "expected_count.dge.txt.gz"
		File stats = "stats.tsv"
	}

	runtime {
		docker: docker
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk " + ceil(size(count_results, "GB")*2 + size(gene_results, "GB")*2 + extra_disk_space)+ " HDD"
		cpu: 1
		preemptible: preemptible
	}
}
