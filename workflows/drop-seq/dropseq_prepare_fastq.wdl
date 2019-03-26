workflow dropseq_prepare_fastq {
	String sample_id
	# A comma-separated list of input FASTQs directories (gs urls)
	String r1
	String r2
	Int preemptible = 2
	Int disk_space
	String? zones = "us-east1-d us-west1-a us-west1-b"
	String output_directory
	String workflow_version

	call PrepareFastq {
		input:
			preemptible=preemptible,
			memory="3750M",
			cpu=1,
			r1=r1,
			r2=r2,
			sample_id=sample_id,
			disk_space=disk_space,
			output_directory=output_directory,
			zones= zones,
			workflow_version=workflow_version
	}

	output {
		String bam=PrepareFastq.bam
		String cellular_tag_summary=PrepareFastq.cellular_tag_summary
      	String molecular_tag_summary=PrepareFastq.molecular_tag_summary
        String adapter_trimming_report = PrepareFastq.adapter_trimming_report
       	String polyA_trimming_report= PrepareFastq.polyA_trimming_report
	}
}


task PrepareFastq {
	String memory
	Int cpu
	String r1
	String r2
	Int disk_space
	String sample_id
	Int preemptible
	String zones
	String output_directory
	String workflow_version

	command {
		set -o pipefail
		set -e

		python <<CODE
		import os
		from subprocess import check_call
		r1_files = '${r1}'.split(',')
		r2_files = '${r2}'.split(',')
		for i in range(len(r1_files)):
			call_args = ['gsutil', '-q', '-m', 'cp', r1_files[i], r2_files[i], '.']
			check_call(call_args)
			call_args = ['mv', os.path.basename(r1_files[i]), 'S_R1_' + str(i+1).zfill(3) + '.fastq.gz']
			check_call(call_args)
			call_args = ['mv', os.path.basename(r2_files[i]), 'S_R2_' + str(i+1).zfill(3) + '.fastq.gz']
			check_call(call_args)

			# rename to fastq files to  S_R1_###.fastq.gz
		CODE

		mkfifo pipe1 pipe2 pipe3 pipe4 pipe5

		java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/picard.jar FastqToSam \
		OUTPUT=pipe1 \
		USE_SEQUENTIAL_FASTQS=true \
		FASTQ=S_R1_001.fastq.gz \
		FASTQ2=S_R2_001.fastq.gz \
		QUALITY_FORMAT=Standard \
		SAMPLE_NAME=${sample_id} \
		SORT_ORDER=queryname | \
		java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
		INPUT=pipe1 \
		OUTPUT=pipe2 \
		SUMMARY=${sample_id}_tagged_cellular_summary.txt \
		BASE_RANGE=1-12 \
		BASE_QUALITY=10 \
		BARCODED_READ=1 \
		DISCARD_READ=false \
		TAG_NAME=XC \
		NUM_BASES_BELOW_QUALITY=1 | \
		java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
		INPUT=pipe2 \
		OUTPUT=pipe3 \
		SUMMARY="${sample_id}_tagged_molecular_summary.txt" \
		BASE_RANGE=13-20 \
		BASE_QUALITY=10 \
		BARCODED_READ=1 \
		DISCARD_READ=true \
		TAG_NAME=XM \
		NUM_BASES_BELOW_QUALITY=1 | \
		java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar FilterBam VALIDATION_STRINGENCY=SILENT \
		INPUT=pipe3 \
		OUTPUT=pipe4 \
		TAG_REJECT=XQ | \
		java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TrimStartingSequence VALIDATION_STRINGENCY=SILENT \
		INPUT=pipe4 \
		OUTPUT=pipe5 \
		OUTPUT_SUMMARY="${sample_id}_adapter_trimming_report.txt" \
		SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
		MISMATCHES=0 \
		NUM_BASES=5 | \
		java -Dsamjdk.compression_level=2 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar PolyATrimmer VALIDATION_STRINGENCY=SILENT \
		INPUT=pipe5 \
		OUTPUT="${sample_id}_aligner_input.bam" \
		OUTPUT_SUMMARY="${sample_id}_polyA_trimming_report.txt" \
		MISMATCHES=0 \
		NUM_BASES=6 \
		NEW=true

		gsutil -q -m cp ${sample_id}_* ${output_directory}/

	}

	output {

		String bam="${output_directory}/${sample_id}_aligner_input.bam"
		String cellular_tag_summary="${output_directory}/${sample_id}_tagged_cellular_summary.txt"
		String molecular_tag_summary="${output_directory}/${sample_id}_tagged_molecular_summary.txt"
		String adapter_trimming_report ="${output_directory}/${sample_id}_adapter_trimming_report.txt"
		String polyA_trimming_report="${output_directory}/${sample_id}_polyA_trimming_report.txt"
	}

	runtime {
		docker: "regevlab/dropseq-${workflow_version}"
		disks: "local-disk ${disk_space} HDD"
		memory :"${memory}"
		cpu:"${cpu}"
		zones: zones
		preemptible: "${preemptible}"
	}
}

