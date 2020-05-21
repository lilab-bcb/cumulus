version 1.0

workflow dropseq_prepare_fastq {
	input {
		String sample_id
		# A comma-separated list of input FASTQs directories (gs urls)
		String r1
		String r2
		Int preemptible = 2
		Int disk_space
		String? zones = "us-east1-d us-west1-a us-west1-b"
		String output_directory
		String drop_seq_tools_version
		Boolean quality_tags
		String umi_base_range
		String cellular_barcode_base_range
		String trim_sequence
		Int trim_num_bases
		String docker_registry
	}
	call FastqToSam {
		input:
			preemptible=preemptible,
			memory="1500M",
			cpu=1,
			r1=r1,
			r2=r2,
			sample_id=sample_id,
			output_directory=output_directory,
			disk_space=disk_space,
			zones= zones,
			drop_seq_tools_version=drop_seq_tools_version,
			docker_registry=docker_registry
	}

	call TagBam {
		input:
			preemptible=preemptible,
			memory="1000M",
			cpu=1,
			unmapped_bam=FastqToSam.bam,
			umi_base_range=umi_base_range,
			cellular_barcode_base_range=cellular_barcode_base_range,
			quality_tags=quality_tags,
			sample_id=sample_id,
			output_directory=output_directory,
			disk_space=disk_space,
			zones= zones,
			drop_seq_tools_version=drop_seq_tools_version,
			docker_registry=docker_registry
	}

	call FilterBam {
		input:
			preemptible=preemptible,
			memory="1000M",
			cpu=1,
			unmapped_bam=TagBam.bam,
			sample_id=sample_id,
			output_directory=output_directory,
			disk_space=disk_space,
			zones= zones,
			drop_seq_tools_version=drop_seq_tools_version,
			docker_registry=docker_registry
	}

	call TrimBam{
		input:
			preemptible=preemptible,
			memory="1000M",
			trim_sequence=trim_sequence,
			trim_num_bases=trim_num_bases,
			cpu=1,
			unmapped_bam=FilterBam.bam,
			sample_id=sample_id,
			output_directory=output_directory,
			disk_space=disk_space,
			zones= zones,
			drop_seq_tools_version=drop_seq_tools_version,
			docker_registry=docker_registry
	}



	output {
		String unmapped_bam=FastqToSam.bam
		String tagged_bam=TagBam.bam
		String filtered_bam=FilterBam.bam
		String trimmed_bam=TrimBam.bam
		String cellular_tag_summary=TagBam.cellular_tag_summary
		String molecular_tag_summary=TagBam.molecular_tag_summary
		String adapter_trimming_report = TrimBam.adapter_trimming_report
		String polyA_trimming_report= TrimBam.polyA_trimming_report
	}
}


task FastqToSam {
	input {
		String memory
		Int cpu
		String r1
		String r2
		Boolean is_gz = sub(r1, "^.+\\.gz$", "GZ") == "GZ"
		String fastq_suffix = if is_gz then ".gz" else ""
		Int disk_space
		String sample_id
		Int preemptible
		String zones
		String output_directory
		String drop_seq_tools_version
		String docker_registry
	}
	command {
		set -o pipefail
		set -e

		

		python <<CODE
		import os
		from subprocess import check_call
		fastq_suffix = '~{fastq_suffix}'
		r1_files = '~{r1}'.split(',')
		r2_files = '~{r2}'.split(',')
		for i in range(len(r1_files)):
			call_args = ['gsutil', '-q', '-m', 'cp', r1_files[i], r2_files[i], '.']
			check_call(call_args)
			call_args = ['mv', os.path.basename(r1_files[i]), 'S_R1_' + str(i+1).zfill(3) + '.fastq' + fastq_suffix]
			check_call(call_args)
			call_args = ['mv', os.path.basename(r2_files[i]), 'S_R2_' + str(i+1).zfill(3) + '.fastq' + fastq_suffix]
			check_call(call_args)
			# rename to fastq files to  S_R1_###.fastq.gz or S_R1_###.fastq
		CODE


		java -Dsamjdk.compression_level=6 -Xmx3000m -jar /software/picard.jar FastqToSam \
		OUTPUT="~{sample_id}_unmapped.bam" \
		USE_SEQUENTIAL_FASTQS=true \
		FASTQ=S_R1_001.fastq~{fastq_suffix} \
		FASTQ2=S_R2_001.fastq~{fastq_suffix} \
		QUALITY_FORMAT=Standard \
		SAMPLE_NAME=~{sample_id} \
		SORT_ORDER=queryname

		gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
	}

	output {
		
		String bam="~{output_directory}/~{sample_id}_unmapped.bam"
	}

	runtime {
		docker: "~{docker_registry}/dropseq:~{drop_seq_tools_version}"
		disks: "local-disk ~{disk_space} HDD"
		memory :"~{memory}"
		cpu:"~{cpu}"
		zones: zones
		preemptible: "~{preemptible}"
	}
}

task TagBam {
	input {
		String memory
		Int cpu
		File unmapped_bam
		Int disk_space
		String sample_id
		Int preemptible
		String zones
		String drop_seq_tools_version
		Boolean quality_tags
		String umi_base_range
		String cellular_barcode_base_range
		String docker_registry
		String output_directory
	}
	command {
		set -o pipefail
		set -e

		
		mkfifo pipe1

		java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
		INPUT=~{unmapped_bam} \
		OUTPUT=pipe1 \
		~{true='BARCODE_QUALITY_TAG=CY' false='' quality_tags} \
		SUMMARY=~{sample_id}_tagged_cellular_summary.txt \
		BASE_RANGE=~{cellular_barcode_base_range} \
		BASE_QUALITY=10 \
		BARCODED_READ=1 \
		DISCARD_READ=false \
		TAG_NAME=XC \
		NUM_BASES_BELOW_QUALITY=1 | \
		java -Dsamjdk.compression_level=6 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
		INPUT=pipe1 \
		OUTPUT="~{sample_id}_tagged.bam" \
		~{true='BARCODE_QUALITY_TAG=UY' false='' quality_tags} \
		SUMMARY="~{sample_id}_tagged_molecular_summary.txt" \
		BASE_RANGE=~{umi_base_range} \
		BASE_QUALITY=10 \
		BARCODED_READ=1 \
		DISCARD_READ=true \
		TAG_NAME=XM \
		NUM_BASES_BELOW_QUALITY=1

		gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
	}

	output {
		
		String bam="~{output_directory}/~{sample_id}_tagged.bam"
		String cellular_tag_summary="~{output_directory}/~{sample_id}_tagged_cellular_summary.txt"
		String molecular_tag_summary="~{output_directory}/~{sample_id}_tagged_molecular_summary.txt"
	}

	runtime {
		docker: "~{docker_registry}/dropseq:~{drop_seq_tools_version}"
		disks: "local-disk ~{disk_space} HDD"
		memory :"~{memory}"
		cpu:"~{cpu}"
		zones: zones
		preemptible: "~{preemptible}"
	}
}

task FilterBam {
	input {
		String memory
		Int cpu
		File unmapped_bam
		Int disk_space
		String sample_id
		Int preemptible
		String zones
		String output_directory
		String drop_seq_tools_version
		String docker_registry
	}
	command {
		set -o pipefail
		set -e

		
		java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar FilterBam VALIDATION_STRINGENCY=SILENT \
				INPUT=~{unmapped_bam} \
				OUTPUT="~{sample_id}_filtered.bam" \
				TAG_REJECT=XQ
		gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
	}

	output {
		
		String bam="~{output_directory}/~{sample_id}_filtered.bam"
	}

	runtime {
		docker: "~{docker_registry}/dropseq:~{drop_seq_tools_version}"
		disks: "local-disk ~{disk_space} HDD"
		memory :"~{memory}"
		cpu:"~{cpu}"
		zones: zones
		preemptible: "~{preemptible}"
	}
}

task TrimBam {
	input {
		String memory
		Int cpu
		File unmapped_bam
		Int disk_space
		String sample_id
		Int preemptible
		String zones
		String output_directory
		String drop_seq_tools_version
		String trim_sequence
		Int trim_num_bases
		String docker_registry
	}
	command {
		set -o pipefail
		set -e

		
		mkfifo pipe1

		java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TrimStartingSequence VALIDATION_STRINGENCY=SILENT \
				INPUT=~{unmapped_bam} \
				OUTPUT=pipe1 \
				OUTPUT_SUMMARY="~{sample_id}_adapter_trimming_report.txt" \
				SEQUENCE=~{trim_sequence} \
				MISMATCHES=0 \
				NUM_BASES=~{trim_num_bases} | \
				java -Dsamjdk.compression_level=2 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar PolyATrimmer VALIDATION_STRINGENCY=SILENT \
				INPUT=pipe1 \
				OUTPUT="~{sample_id}_trimmed.bam" \
				OUTPUT_SUMMARY="~{sample_id}_polyA_trimming_report.txt" \
				MISMATCHES=0 \
				NUM_BASES=6 \
				NEW=true

		gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
	}

	output {
		
		String bam="~{output_directory}/~{sample_id}_trimmed.bam"
		String adapter_trimming_report ="~{output_directory}/~{sample_id}_adapter_trimming_report.txt"
		String polyA_trimming_report="~{output_directory}/~{sample_id}_polyA_trimming_report.txt"
	}

	runtime {
		docker: "~{docker_registry}/dropseq:~{drop_seq_tools_version}"
		disks: "local-disk ~{disk_space} HDD"
		memory :"~{memory}"
		cpu:"~{cpu}"
		zones: zones
		preemptible: "~{preemptible}"
	}
}


