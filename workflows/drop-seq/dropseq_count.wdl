workflow dropseq_count {
	String sample_id
	File input_bam
	Int? force_cells
	Int preemptible = 2
	String? zones = "us-east1-d us-west1-a us-west1-b"
	String output_directory
	String drop_seq_tools_version
	File? cellular_barcode_whitelist
	Float? disk_multiplier = 3.5

	call DigitalExpressionPrep {
		input:
			preemptible=preemptible,
			output_directory=output_directory,
			input_bam=input_bam,
			sample_id=sample_id,
			cellular_barcode_whitelist=cellular_barcode_whitelist,
			memory="3750M",
			disk_multiplier=disk_multiplier,
			zones=zones,
			drop_seq_tools_version=drop_seq_tools_version
	}

	call CollectCellBarcodes {
		input:
			preemptible=preemptible,
			output_directory=output_directory,
			histogram=DigitalExpressionPrep.histogram,
			ncells=if defined(force_cells) then force_cells	 else read_int(DigitalExpressionPrep.ncells),
			sample_id=sample_id,
			zones=zones,
			drop_seq_tools_version=drop_seq_tools_version
	}

	call DigitalExpression {
		input:
			preemptible=preemptible,
			output_directory=output_directory,
			input_bam=DigitalExpressionPrep.bam,
			sample_id=sample_id,
			barcodes=CollectCellBarcodes.cell_barcodes,
			memory="3750M",
			zones=zones,
			drop_seq_tools_version=drop_seq_tools_version
	}
	output {
		String bead_synthesis_stats = DigitalExpressionPrep.bead_synthesis_stats
		String bead_synthesis_summary = DigitalExpressionPrep.bead_synthesis_summary
		String bead_synthesis_report = DigitalExpressionPrep.bead_synthesis_report
		String bead_substitution_report =DigitalExpressionPrep.bead_substitution_report
		String histogram=DigitalExpressionPrep.histogram
		String ncells=DigitalExpressionPrep.ncells
		String cumplot=DigitalExpressionPrep.cumplot
		String reads_plot=DigitalExpressionPrep.reads_plot

		String cell_barcodes = CollectCellBarcodes.cell_barcodes

		String dge=DigitalExpression.dge
		String dge_summary = DigitalExpression.dge_summary
		String dge_reads=DigitalExpression.dge_reads
		String dge_summary_reads = DigitalExpression.dge_summary_reads
	}
	
}

task CollectCellBarcodes {
	File histogram
	Int ncells
	String sample_id
	Int preemptible
	String zones
	String output_directory
	String drop_seq_tools_version

	command {
		set -e

		collect_cell_barcodes.R \
		${histogram} \
		${ncells} \
		${sample_id}_barcodes_use.txt

		gsutil -m -q cp "${sample_id}_barcodes_use.txt" ${output_directory}/
	}

	output {

		String cell_barcodes="${output_directory}/${sample_id}_barcodes_use.txt"
	}
	runtime {
		docker: "regevlab/dropseq-${drop_seq_tools_version}"
		zones: zones
		disks: "local-disk " + sub(((size(histogram,"GB")+1)*2),"\\..*","") + " HDD"
		preemptible: "${preemptible}"
	}
}



task DigitalExpression {
	String memory
	File input_bam
	File barcodes
	String sample_id
	Int preemptible
	String zones
	String output_directory
	String drop_seq_tools_version


	command {
		set -e

		# CODING+UTR regions on the SENSE strand by default
		java -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar DigitalExpression VALIDATION_STRINGENCY=SILENT \
		I=${input_bam} \
		O="${sample_id}_dge.txt.gz" \
		SUMMARY="${sample_id}_dge.summary.txt" \
		CELL_BC_FILE=${barcodes}

		java -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar DigitalExpression VALIDATION_STRINGENCY=SILENT \
		I=${input_bam} \
		O="${sample_id}_dge_reads.txt.gz" \
		SUMMARY="${sample_id}_dge_reads.summary.txt" \
		CELL_BC_FILE=${barcodes} \
		OUTPUT_READS_INSTEAD=true

		gsutil -q -m cp ${sample_id}_* ${output_directory}/
	}

	output {

		String dge="${output_directory}/${sample_id}_dge.txt.gz"
		String dge_summary = "${output_directory}/${sample_id}_dge.summary.txt"
		String dge_reads="${output_directory}/${sample_id}_dge_reads.txt.gz"
		String dge_summary_reads = "${output_directory}/${sample_id}_dge_reads.summary.txt"
	}

	runtime {
		docker: "regevlab/dropseq-${drop_seq_tools_version}"
		disks: "local-disk " + sub(((size(input_bam,"GB")+1)*3.25),"\\..*","") + " HDD"
		memory :"${memory}"
		zones: zones
		preemptible: "${preemptible}"
		cpu:1
	}
}

task DigitalExpressionPrep {
	String memory
	File input_bam
	String sample_id
	Int preemptible
	String zones
	String output_directory
	String drop_seq_tools_version
	File? cellular_barcode_whitelist
	Float disk_multiplier
	command {
		set -e

		java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar DetectBeadSynthesisErrors VALIDATION_STRINGENCY=SILENT \
		INPUT=${input_bam} \
		MIN_UMIS_PER_CELL=20 \
		OUTPUT_STATS=${sample_id}_synthesis_error_stats.txt \
		SUMMARY=${sample_id}_synthesis_error_summary.txt \
		REPORT=${sample_id}_synthesis_error_report.txt \
		OUTPUT=temp.bam

		java -Dsamjdk.compression_level=2 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar DetectBeadSubstitutionErrors VALIDATION_STRINGENCY=SILENT \
		INPUT=temp.bam \
		OUTPUT="${sample_id}_aligned_tagged_repaired.bam" \
		MIN_UMIS_PER_CELL=20 \
		OUTPUT_REPORT="${sample_id}_substitution_error_report.txt"

		java -Dsamjdk.compression_level=2 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar BamTagHistogram VALIDATION_STRINGENCY=SILENT \
		I="${sample_id}_aligned_tagged_repaired.bam" \
		O=${sample_id}_tag.txt.gz \
		TAG=XC

		python <<CODE
		from subprocess import check_call
		whitelist = '${cellular_barcode_whitelist}'
		if whitelist !='':
			check_call(['filter_histogram.py', '--histogram', '"${sample_id}_tag.txt.gz"', '--whitelist', '${cellular_barcode_whitelist}', '--output', '"${sample_id}_tag.txt.gz"'])
		CODE

		DropSeqCumuPlot.R \
		--collapsed "${sample_id}_tag.txt.gz" \
		--counts "${sample_id}_tag.txt.gz" \
		--num_cells ${sample_id}_ncells.txt \
		--cummulative_plot ${sample_id}_cumplot.pdf \
		--reads_plot ${sample_id}_NreadsHiToLo.pdf

		gsutil -q -m cp ${sample_id}_* ${output_directory}/
	}

	output {

		String bam="${output_directory}/${sample_id}_aligned_tagged_repaired.bam"
		String bead_synthesis_stats = "${output_directory}/${sample_id}_synthesis_error_stats.txt"
		String bead_synthesis_summary = "${output_directory}/${sample_id}_synthesis_error_summary.txt"
		String bead_synthesis_report = "${output_directory}/${sample_id}_synthesis_error_report.txt"
		String bead_substitution_report = "${output_directory}/${sample_id}_substitution_error_report.txt"
		String histogram="${output_directory}/${sample_id}_tag.txt.gz"
		String? unfiltered_histogram="${output_directory}/${sample_id}_unfiltered_tag.txt.gz"
		String ncells="${output_directory}/${sample_id}_ncells.txt"
		String cumplot="${output_directory}/${sample_id}_cumplot.pdf"
		String reads_plot="${output_directory}/${sample_id}_NreadsHiToLo.pdf"
	}

	runtime {
		docker: "regevlab/dropseq-${drop_seq_tools_version}"
		disks: "local-disk " + ceil(disk_multiplier * size(input_bam, "GB") + 20)+ " HDD"
		memory :"${memory}"
		zones: zones
		preemptible: "${preemptible}"
		cpu:1
	}
}







