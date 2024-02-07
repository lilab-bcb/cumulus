version 1.0


workflow dropseq_count_multi_species {
  input {
    String dge_memory = "3750M"
    String drop_seq_tools_version
    String sample_id
    File barcodes
    Int preemptible = 2
    String zones = "us-east1-d us-west1-a us-west1-b"
    String docker_registry
    File input_bam
    String species
    String output_directory
  }
  call DigitalExpression {
    input:
      memory = dge_memory,
      input_bam = FilterBamBySpecies.bam,
      barcodes = barcodes,
      sample_id = sample_id + "_" + species,
      preemptible = preemptible,
      zones = zones,
      output_directory = output_directory,
      drop_seq_tools_version = drop_seq_tools_version,
      docker_registry = docker_registry
  }
  call FilterBamBySpecies {
    input:
      memory = "1000M",
      cpu = 1,
      input_bam = input_bam,
      species = species,
      preemptible = preemptible,
      zones = zones,
      drop_seq_tools_version = drop_seq_tools_version,
      docker_registry = docker_registry
  }

  output {
    String dge = DigitalExpression.dge
    String dge_reads = DigitalExpression.dge_reads
    String dge_summary = DigitalExpression.dge_summary
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


  output {
    String dge = "${output_directory}/${sample_id}_dge.txt.gz"
    String dge_summary = "${output_directory}/${sample_id}_dge.summary.txt"
    String dge_reads = "${output_directory}/${sample_id}_dge_reads.txt.gz"
    String dge_summary_reads = "${output_directory}/${sample_id}_dge_reads.summary.txt"
  }
  command <<<

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
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil(((size(input_bam, "GB") + 1) * 3.25)) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "${memory}"
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


  output {
    File bam = "${species}.bam"
  }
  command <<<


      set -e

      java -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar FilterBam VALIDATION_STRINGENCY=SILENT \
              INPUT=~{input_bam} \
              OUTPUT="~{species}.bam" \
              REF_SOFT_MATCHED_RETAINED=~{species}

   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil(size(input_bam, "GB") * 3 + 5) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
  }

}

