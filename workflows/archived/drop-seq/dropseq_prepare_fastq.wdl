version 1.0


workflow dropseq_prepare_fastq {
  input {
    Int disk_space
    String fastq_to_sam_memory
    String trim_bam_memory
    Int preemptible = 2
    String umi_base_range
    String docker_registry
    String r1
    String zones = "us-east1-d us-west1-a us-west1-b"
    String cellular_barcode_base_range
    String sample_id
    String drop_seq_tools_version
    String r2
    Int trim_num_bases
    String trim_sequence
    Boolean quality_tags
    String output_directory
  }
  call FilterBam {
    input:
      memory = "1000M",
      cpu = 1,
      unmapped_bam = TagBam.bam,
      disk_space = disk_space,
      sample_id = sample_id,
      preemptible = preemptible,
      zones = zones,
      output_directory = output_directory,
      drop_seq_tools_version = drop_seq_tools_version,
      docker_registry = docker_registry
  }
  call FastqToSam {
    input:
      memory = fastq_to_sam_memory,
      cpu = 1,
      r1 = r1,
      r2 = r2,
      disk_space = disk_space,
      sample_id = sample_id,
      preemptible = preemptible,
      zones = zones,
      output_directory = output_directory,
      drop_seq_tools_version = drop_seq_tools_version,
      docker_registry = docker_registry
  }
  call TrimBam {
    input:
      memory = trim_bam_memory,
      cpu = 1,
      unmapped_bam = FilterBam.bam,
      disk_space = disk_space,
      sample_id = sample_id,
      preemptible = preemptible,
      zones = zones,
      output_directory = output_directory,
      drop_seq_tools_version = drop_seq_tools_version,
      trim_sequence = trim_sequence,
      trim_num_bases = trim_num_bases,
      docker_registry = docker_registry
  }
  call TagBam {
    input:
      memory = "1000M",
      cpu = 1,
      unmapped_bam = FastqToSam.bam,
      disk_space = disk_space,
      sample_id = sample_id,
      preemptible = preemptible,
      zones = zones,
      drop_seq_tools_version = drop_seq_tools_version,
      quality_tags = quality_tags,
      umi_base_range = umi_base_range,
      cellular_barcode_base_range = cellular_barcode_base_range,
      docker_registry = docker_registry,
      output_directory = output_directory
  }

  output {
    String unmapped_bam = FastqToSam.bam
    String cellular_tag_summary = TagBam.cellular_tag_summary
    String polyA_trimming_report = TrimBam.polyA_trimming_report
    String filtered_bam = FilterBam.bam
    String trimmed_bam = TrimBam.bam
    String molecular_tag_summary = TagBam.molecular_tag_summary
    String tagged_bam = TagBam.bam
    String adapter_trimming_report = TrimBam.adapter_trimming_report
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


  output {
    String bam = "${output_directory}/${sample_id}_filtered.bam"
  }
  command <<<

      set -o pipefail
      set -e

        
      java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar FilterBam VALIDATION_STRINGENCY=SILENT \
              INPUT=~{unmapped_bam} \
              OUTPUT="~{sample_id}_filtered.bam" \
              TAG_REJECT=XQ
        gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk ${disk_space} HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
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


  output {
    String bam = "${output_directory}/${sample_id}_tagged.bam"
    String cellular_tag_summary = "${output_directory}/${sample_id}_tagged_cellular_summary.txt"
    String molecular_tag_summary = "${output_directory}/${sample_id}_tagged_molecular_summary.txt"
  }
  command <<<

      set -o pipefail
      set -e

        
      mkfifo pipe1

      java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
      INPUT=~{unmapped_bam} \
      OUTPUT=pipe1 \
      ~{true="BARCODE_QUALITY_TAG=CY" false=""  quality_tags} \
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
      ~{true="BARCODE_QUALITY_TAG=UY" false=""  quality_tags} \
      SUMMARY="~{sample_id}_tagged_molecular_summary.txt" \
      BASE_RANGE=~{umi_base_range} \
      BASE_QUALITY=10 \
      BARCODED_READ=1 \
      DISCARD_READ=true \
      TAG_NAME=XM \
      NUM_BASES_BELOW_QUALITY=1

      gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk ${disk_space} HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
  }

}
task FastqToSam {
  input {
    String memory
    Int cpu
    String r1
    String r2
    Int disk_space
    String sample_id
    Int preemptible
    String zones
    String output_directory
    String drop_seq_tools_version
    String docker_registry
  }
    Boolean is_gz = sub(r1, "^.+\.gz$", "GZ") == "GZ"
    String fastq_suffix = if is_gz then ".gz" else ""



  output {
    String bam = "${output_directory}/${sample_id}_unmapped.bam"
  }
  command <<<

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

      command_mem_mb=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2) }')
      command_mem_mb=$(($command_mem_mb-200))
      command_mem_mb=$command_mem_mb"m"

      java -Dsamjdk.compression_level=6 -Xmx$command_mem_mb -jar /software/picard.jar FastqToSam \
      OUTPUT="~{sample_id}_unmapped.bam" \
      USE_SEQUENTIAL_FASTQS=true \
      FASTQ=S_R1_001.fastq~{fastq_suffix} \
      FASTQ2=S_R2_001.fastq~{fastq_suffix} \
      QUALITY_FORMAT=Standard \
      SAMPLE_NAME=~{sample_id} \
      SORT_ORDER=queryname

      gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk ${disk_space} HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
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

  output {

    String bam = "${output_directory}/${sample_id}_trimmed.bam"
    String adapter_trimming_report = "${output_directory}/${sample_id}_adapter_trimming_report.txt"
    String polyA_trimming_report = "${output_directory}/${sample_id}_polyA_trimming_report.txt"
  }
  command <<<

        set -o pipefail
        set -e
        
        mkfifo pipe1
        command_mem_mb=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2) }')
        command_mem_mb=$(($command_mem_mb-200))
        command_mem_mb=$command_mem_mb"m"

        java -Dsamjdk.compression_level=1 -Xmx$command_mem_mb -jar /software/Drop-seq_tools/jar/dropseq.jar TrimStartingSequence VALIDATION_STRINGENCY=SILENT \
              INPUT=~{unmapped_bam} \
              OUTPUT=pipe1 \
              OUTPUT_SUMMARY="~{sample_id}_adapter_trimming_report.txt" \
              SEQUENCE=~{trim_sequence} \
              MISMATCHES=0 \
              NUM_BASES=~{trim_num_bases} | \
              java -Dsamjdk.compression_level=2 -Xmx$command_mem_mb -jar /software/Drop-seq_tools/jar/dropseq.jar PolyATrimmer VALIDATION_STRINGENCY=SILENT \
              INPUT=pipe1 \
              OUTPUT="~{sample_id}_trimmed.bam" \
              OUTPUT_SUMMARY="~{sample_id}_polyA_trimming_report.txt" \
              MISMATCHES=0 \
              NUM_BASES=6 \
              NEW=true

        gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
    
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk ${disk_space} HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
  }

}

