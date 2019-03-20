workflow dropseq_qc {
   	String sample_id
    File input_bam
    File cell_barcodes
	File refflat
	String? zones = "us-east1-d us-west1-a us-west1-b"
	Int preemptible = 2
	String output_directory
	String workflow_version

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
			workflow_version=workflow_version
	}


	output {
		String sc_rnaseq_metrics_report=SingleCellRnaSeqMetricsCollector.sc_rnaseq_metrics_report
	}

}


task SingleCellRnaSeqMetricsCollector {
	 File input_bam
     File refflat
     File cell_barcodes
     String sample_id
     String output_directory
     String zones
     Int preemptible
     String workflow_version
     String memory



	command {
    	set -e
		monitor_script.sh > monitoring.log &

		java -Xmx3000m -jar /software/Drop-seq_tools-2.1.0/jar/dropseq.jar SingleCellRnaSeqMetricsCollector VALIDATION_STRINGENCY=SILENT \
		INPUT=${input_bam} \
		OUTPUT=${sample_id}_rnaseq_metrics.txt \
		ANNOTATIONS_FILE=${refflat} \
		CELL_BC_FILE=${cell_barcodes}

		gsutil -q -m cp *.txt ${output_directory}/
	}

	output {
		File monitoringLog = "monitoring.log"
		File sc_rnaseq_metrics_report = "${sample_id}_rnaseq_metrics.txt"
	}

	runtime {
		docker: "regevlab/dropseq-${workflow_version}"
		disks: "local-disk " + ceil(2 + size(input_bam,"GB") + size(cell_barcodes, "GB") + size(refflat,"GB")) + " HDD"
		memory :"${memory}"
		preemptible: "${preemptible}"
		zones: zones
		cpu:1
	}
}







