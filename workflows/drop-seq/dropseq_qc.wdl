workflow dropseq_qc {
   	String sample_id
    File input_bam
    File? cell_barcodes
	File refflat
	String? zones = "us-east1-d us-west1-a us-west1-b"
	Int preemptible = 2
	String output_directory
	String drop_seq_tools_version
	String docker_registry

	call SingleCellRnaSeqMetricsCollector {
		input:
			memory="3750M",
			input_bam=input_bam,
			sample_id=sample_id,
			cell_barcodes=cell_barcodes,
			preemptible=preemptible,
			refflat=refflat,
			zones=zones,
			output_directory=output_directory,
			drop_seq_tools_version=drop_seq_tools_version,
			docker_registry=docker_registry
	}


	output {
		String sc_rnaseq_metrics_report=SingleCellRnaSeqMetricsCollector.sc_rnaseq_metrics_report
	}

}


task SingleCellRnaSeqMetricsCollector {
    File input_bam
    File refflat
    File? cell_barcodes
    String sample_id
    String output_directory
    String zones
    Int preemptible
    String drop_seq_tools_version
    String memory
    String docker_registry

	command {
    	set -e

		java -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar SingleCellRnaSeqMetricsCollector VALIDATION_STRINGENCY=SILENT \
		INPUT=${input_bam} \
		OUTPUT=${sample_id}_rnaseq_metrics.txt \
		ANNOTATIONS_FILE=${refflat} \
		${"CELL_BC_FILE=" + cell_barcodes}

		gsutil -q -m cp *.txt ${output_directory}/
	}

	output {

		File sc_rnaseq_metrics_report = "${sample_id}_rnaseq_metrics.txt"
	}

	runtime {
		docker: "${docker_registry}dropseq:${drop_seq_tools_version}"
		disks: "local-disk " + ceil(20 + 3.25*size(input_bam,"GB") + size(cell_barcodes, "GB") + size(refflat,"GB")) + " HDD"
		memory :"${memory}"
		preemptible: "${preemptible}"
		zones: zones
		cpu:1
	}
}







