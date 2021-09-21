version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropseq_count_multi_species/versions/2/plain-WDL/descriptor" as dropseq_count_multi_species_wdl
#import "dropseq_count_multi_species.wdl" as dropseq_count_multi_species_wdl

workflow dropseq_count {
  input {
    Float disk_multiplier = 3.5
    Int preemptible = 2
    File? cellular_barcode_whitelist
    String dge_prep_memory = "3750M"
    String sample_id
    String dge_memory = "3750M"
    String output_directory
    String docker_registry
    String zones = "us-east1-d us-west1-a us-west1-b"
    Int? force_cells
    String drop_seq_tools_version
    Array[String] species
    File input_bam
  }
  if (!defined(species) || length(species) <= 1) {
    call DigitalExpression {
      input:
        memory = dge_memory,
        input_bam = DigitalExpressionPrep.bam,
        barcodes = CollectCellBarcodes.cell_barcodes,
        sample_id = sample_id,
        preemptible = preemptible,
        zones = zones,
        output_directory = output_directory,
        drop_seq_tools_version = drop_seq_tools_version,
        docker_registry = docker_registry
    }
  }
  if (defined(species) && length(species) > 1) {
    scatter (s in species) {
      call dropseq_count_multi_species_wdl.dropseq_count_multi_species as dge_species {
        input:
          dge_memory = dge_memory,
          sample_id = sample_id,
          preemptible = preemptible,
          zones = zones,
          output_directory = output_directory,
          docker_registry = docker_registry,
          drop_seq_tools_version = drop_seq_tools_version,
          barcodes = CollectCellBarcodes.cell_barcodes,
          species = s,
          input_bam = input_bam
      }
    }
  }
  call DigitalExpressionPrep {
    input:
      memory = dge_prep_memory,
      input_bam = input_bam,
      sample_id = sample_id,
      preemptible = preemptible,
      zones = zones,
      output_directory = output_directory,
      drop_seq_tools_version = drop_seq_tools_version,
      cellular_barcode_whitelist = cellular_barcode_whitelist,
      disk_multiplier = disk_multiplier,
      docker_registry = docker_registry
  }
  call CollectCellBarcodes {
    input:
      histogram = DigitalExpressionPrep.histogram,
      ncells = if defined(force_cells) then force_cells else read_int(DigitalExpressionPrep.ncells),
      sample_id = sample_id,
      preemptible = preemptible,
      zones = zones,
      output_directory = output_directory,
      drop_seq_tools_version = drop_seq_tools_version,
      docker_registry = docker_registry
  }

  output {
    String reads_plot = DigitalExpressionPrep.reads_plot
    String bead_synthesis_report = DigitalExpressionPrep.bead_synthesis_report
    String ncells = DigitalExpressionPrep.ncells
    String cumplot = DigitalExpressionPrep.cumplot
    String? dge_summary = DigitalExpression.dge_summary
    String? dge_reads = DigitalExpression.dge_reads
    String histogram = DigitalExpressionPrep.histogram
    String bead_synthesis_stats = DigitalExpressionPrep.bead_synthesis_stats
    String? dge_summary_reads = DigitalExpression.dge_summary_reads
    String cell_barcodes = CollectCellBarcodes.cell_barcodes
    Array[String]? dge_summary_reads_multi_species = dge_species.dge_summary_reads
    String bead_substitution_report = DigitalExpressionPrep.bead_substitution_report
    Array[String]? dge_reads_multi_species = dge_species.dge_reads
    Array[String]? dge_summary_multi_species = dge_species.dge_summary
    String? dge = DigitalExpression.dge
    Array[String]? dge_multi_species = dge_species.dge
    String bead_synthesis_summary = DigitalExpressionPrep.bead_synthesis_summary
  }
}
task CollectCellBarcodes {
  input {
    File histogram
    Int? ncells
    String sample_id
    Int preemptible
    String zones
    String output_directory
    String drop_seq_tools_version
    String docker_registry
  }


  output {
    String cell_barcodes = "${output_directory}/${sample_id}_barcodes_use.txt"
  }
  command <<<

      set -e

      collect_cell_barcodes.R \
      ~{histogram} \
      ~{ncells} \
      ~{sample_id}_barcodes_use.txt

      gsutil -m -q cp "~{sample_id}_barcodes_use.txt" ~{output_directory}/
   
  >>>
  runtime {
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    zones: zones
    disks: "local-disk " + ceil(((size(histogram, "GB") + 1) * 2)) + " HDD"
    preemptible: "${preemptible}"
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
task DigitalExpressionPrep {
  input {
    String memory
    File input_bam
    String sample_id
    Int preemptible
    String zones
    String output_directory
    String drop_seq_tools_version
    File? cellular_barcode_whitelist
    Float disk_multiplier
    String docker_registry
  }


  output {
    String bam = "${output_directory}/${sample_id}_aligned_tagged_repaired.bam"
    String bead_synthesis_stats = "${output_directory}/${sample_id}_synthesis_error_stats.txt"
    String bead_synthesis_summary = "${output_directory}/${sample_id}_synthesis_error_summary.txt"
    String bead_synthesis_report = "${output_directory}/${sample_id}_synthesis_error_report.txt"
    String bead_substitution_report = "${output_directory}/${sample_id}_substitution_error_report.txt"
    String histogram = "${output_directory}/${sample_id}_tag.txt.gz"
    String? unfiltered_histogram = "${output_directory}/${sample_id}_unfiltered_tag.txt.gz"
    String ncells = "${output_directory}/${sample_id}_ncells.txt"
    String cumplot = "${output_directory}/${sample_id}_cumplot.pdf"
    String reads_plot = "${output_directory}/${sample_id}_NreadsHiToLo.pdf"
  }
  command <<<

      set -e

      java -Dsamjdk.compression_level=1 -Xmx~{memory} -jar /software/Drop-seq_tools/jar/dropseq.jar DetectBeadSynthesisErrors VALIDATION_STRINGENCY=SILENT \
      INPUT=~{input_bam} \
      MIN_UMIS_PER_CELL=20 \
      OUTPUT_STATS=~{sample_id}_synthesis_error_stats.txt \
      SUMMARY=~{sample_id}_synthesis_error_summary.txt \
      REPORT=~{sample_id}_synthesis_error_report.txt \
      OUTPUT=temp.bam

      java -Dsamjdk.compression_level=2 -Xmx~{memory} -jar /software/Drop-seq_tools/jar/dropseq.jar DetectBeadSubstitutionErrors VALIDATION_STRINGENCY=SILENT \
      INPUT=temp.bam \
      OUTPUT="~{sample_id}_aligned_tagged_repaired.bam" \
      MIN_UMIS_PER_CELL=20 \
      OUTPUT_REPORT="~{sample_id}_substitution_error_report.txt"

      java -Dsamjdk.compression_level=2 -Xmx~{memory} -jar /software/Drop-seq_tools/jar/dropseq.jar BamTagHistogram VALIDATION_STRINGENCY=SILENT \
      I="~{sample_id}_aligned_tagged_repaired.bam" \
      O=~{sample_id}_tag.txt.gz \
      TAG=XC

      python <<CODE
      from subprocess import check_call
      whitelist = '~{cellular_barcode_whitelist}'
      if whitelist !='':
         check_call(['filter_histogram.py', '--histogram', '~{sample_id}_tag.txt.gz', '--whitelist', '~{cellular_barcode_whitelist}', '--output', '~{sample_id}_tag.txt.gz'])
      CODE

      DropSeqCumuPlot.R \
      --collapsed "~{sample_id}_tag.txt.gz" \
      --counts "~{sample_id}_tag.txt.gz" \
      --num_cells ~{sample_id}_ncells.txt \
      --cummulative_plot ~{sample_id}_cumplot.pdf \
      --reads_plot ~{sample_id}_NreadsHiToLo.pdf

      gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil(disk_multiplier * size(input_bam, "GB") + 20) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "${memory}"
  }

}

