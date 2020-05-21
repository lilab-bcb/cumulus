version 1.0

workflow dropseq_count_multi_species {
	input {
		String sample_id
		File input_bam
		String species
		String output_directory
		File barcodes

		String? zones = "us-east1-d us-west1-a us-west1-b"
		String? dge_memory = "3750M"
		Int? preemptible = 2
		String docker_registry
		String drop_seq_tools_version
	}
	call FilterBamBySpecies {
		input:
			preemptible=preemptible,
			input_bam=input_bam,
			species=species,
			memory="1000M",
			cpu=1,
			zones=zones,
			drop_seq_tools_version=drop_seq_tools_version,
			docker_registry=docker_registry
	}

	call DigitalExpression {
		input:
			preemptible=preemptible,
			output_directory=output_directory,
			input_bam=FilterBamBySpecies.bam,
			sample_id=sample_id + '_' + species,
			barcodes=barcodes,
			memory=dge_memory,
			zones=zones,
			drop_seq_tools_version=drop_seq_tools_version,
			docker_registry=docker_registry
	}

	output {
		String dge=DigitalExpression.dge
		String dge_summary = DigitalExpression.dge_summary
		String dge_reads=DigitalExpression.dge_reads
		String dge_summary_reads = DigitalExpression.dge_summary_reads
	}
}


task DigitalExpression {
	input {
		String memory
		File input_bam
		File barcodes
		String sample_id
		Int preemptible
		String zones
		String output_directory
		String drop_seq_tools_version
		String docker_registry
	}
	command {
		set -e

		# CODING+UTR regions on the SENSE strand by default
		java -Xmx~{memory} -jar /software/Drop-seq_tools/jar/dropseq.jar DigitalExpression VALIDATION_STRINGENCY=SILENT \
		I=~{input_bam} \
		O="~{sample_id}_dge.txt.gz" \
		SUMMARY="~{sample_id}_dge.summary.txt" \
		CELL_BC_FILE=~{barcodes}

		java -Xmx~{memory} -jar /software/Drop-seq_tools/jar/dropseq.jar DigitalExpression VALIDATION_STRINGENCY=SILENT \
		I=~{input_bam} \
		O="~{sample_id}_dge_reads.txt.gz" \
		SUMMARY="~{sample_id}_dge_reads.summary.txt" \
		CELL_BC_FILE=~{barcodes} \
		OUTPUT_READS_INSTEAD=true

		gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
	}

	output {
		String dge="~{output_directory}/~{sample_id}_dge.txt.gz"
		String dge_summary = "~{output_directory}/~{sample_id}_dge.summary.txt"
		String dge_reads="~{output_directory}/~{sample_id}_dge_reads.txt.gz"
		String dge_summary_reads = "~{output_directory}/~{sample_id}_dge_reads.summary.txt"
	}

	runtime {
		docker: "~{docker_registry}/dropseq:~{drop_seq_tools_version}"
		disks: "local-disk " + sub(((size(input_bam,"GB")+1)*3.25),"\\..*","") + " HDD"
		memory :"~{memory}"
		zones: zones
		preemptible: "~{preemptible}"
		cpu:1
	}
}


task FilterBamBySpecies {
	input {
		String memory
		Int cpu
		File input_bam
		String species
		Int preemptible
		String zones
		String drop_seq_tools_version
		String docker_registry
	}
	command {

		set -e

		java -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar FilterBam VALIDATION_STRINGENCY=SILENT \
				INPUT=~{input_bam} \
				OUTPUT="~{species}.bam" \
				REF_SOFT_MATCHED_RETAINED=~{species}

	}

	output {
		File bam="~{species}.bam"
	}

	runtime {
		docker: "~{docker_registry}/dropseq:~{drop_seq_tools_version}"
		disks: "local-disk " + ceil(size(input_bam, "GB")*3  + 5)+ " HDD"
		memory :"~{memory}"
		cpu:"~{cpu}"
		zones: zones
		preemptible: "~{preemptible}"
	}
}








