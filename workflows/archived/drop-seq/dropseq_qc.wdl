version 1.0


workflow dropseq_qc {
  input {
    String drop_seq_tools_version
    Int preemptible = 2
    String zones = "us-east1-d us-west1-a us-west1-b"
    File? cell_barcodes
    File refflat
    String sample_id
    String output_directory
    File input_bam
    String docker_registry
  }
  call SingleCellRnaSeqMetricsCollector {
    input:
      input_bam = input_bam,
      refflat = refflat,
      cell_barcodes = cell_barcodes,
      sample_id = sample_id,
      output_directory = output_directory,
      zones = zones,
      preemptible = preemptible,
      drop_seq_tools_version = drop_seq_tools_version,
      memory = "3750M",
      docker_registry = docker_registry
  }

  output {
    String sc_rnaseq_metrics_report = SingleCellRnaSeqMetricsCollector.sc_rnaseq_metrics_report
  }
}
task SingleCellRnaSeqMetricsCollector {
  input {
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
  }


  output {
    File sc_rnaseq_metrics_report = "${sample_id}_rnaseq_metrics.txt"
  }
  command <<<

       set -e

      java -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar SingleCellRnaSeqMetricsCollector VALIDATION_STRINGENCY=SILENT \
      INPUT=~{input_bam} \
      OUTPUT=~{sample_id}_rnaseq_metrics.txt \
      ANNOTATIONS_FILE=~{refflat} \
      ~{"CELL_BC_FILE=" + cell_barcodes}

      gsutil -q -m cp *.txt ~{output_directory}/
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil(20 + 3.25 * size(input_bam, "GB") + size(cell_barcodes, "GB") + size(refflat, "GB")) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "${memory}"
  }

}

