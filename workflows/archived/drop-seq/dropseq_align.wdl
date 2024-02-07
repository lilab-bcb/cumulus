version 1.0


workflow dropseq_align {
  input {
    String zones = "us-east1-d us-west1-a us-west1-b"
    Float add_bam_tags_disk_space_multiplier = 25
    String merge_bam_alignment_memory = "13G"
    String drop_seq_tools_version
    String star_memory = "57.6G"
    Int sort_bam_max_records_in_ram = 2000000
    Int star_cpus = 64
    File genome_dict
    String docker_registry
    File gene_intervals
    String sample_id
    File input_bam
    Int preemptible = 2
    String? star_flags
    File refflat
    Float star_extra_disk_space = 2
    File star_genome_file
    String output_directory
    Float star_disk_space_multiplier = 4
    File genome_fasta
  }
  call STAR {
    input:
      memory = star_memory,
      flags = star_flags,
      input_bam = input_bam,
      cpu = star_cpus,
      sample_id = sample_id,
      genome_tar = star_genome_file,
      preemptible = preemptible,
      zones = zones,
      output_directory = output_directory,
      drop_seq_tools_version = drop_seq_tools_version,
      docker_registry = docker_registry,
      disk_space_multiplier = star_disk_space_multiplier,
      extra_disk_space = star_extra_disk_space
  }
  call AddTags {
    input:
      memory = merge_bam_alignment_memory,
      aligned_bam = STAR.bam,
      unaligned_bam = input_bam,
      sample_id = sample_id,
      preemptible = preemptible,
      genome_fasta = genome_fasta,
      genome_dict = genome_dict,
      zones = zones,
      refflat = refflat,
      gene_intervals = gene_intervals,
      output_directory = output_directory,
      drop_seq_tools_version = drop_seq_tools_version,
      disk_space_multiplier = add_bam_tags_disk_space_multiplier,
      cpu = 1,
      sort_bam_max_records_in_ram = sort_bam_max_records_in_ram,
      docker_registry = docker_registry
  }

  output {
    String output_sample_id = AddTags.output_sample_id
    String aligned_bam = STAR.bam
    String star_log_final = STAR.star_log_final
    String aligned_tagged_bam = AddTags.bam
  }
}
task STAR {
  input {
    String memory
    String? flags
    File input_bam
    Int cpu
    String sample_id
    File genome_tar
    Int preemptible
    String zones
    String output_directory
    String drop_seq_tools_version
    String docker_registry
    Float disk_space_multiplier
    Float extra_disk_space
  }


  output {
    String bam = "${output_directory}/${sample_id}_Aligned.out.bam"
    String star_log_final = "${output_directory}/${sample_id}_Log.final.out"
  }
  command <<<

      set -o pipefail
      set -e

      mkdir -p genome_dir
      tar xf ~{genome_tar} -C genome_dir --strip-components 1

      java -Xmx5000m -jar /software/picard.jar SamToFastq \
      VALIDATION_STRINGENCY=SILENT \
      INPUT=~{input_bam} \
      FASTQ=/dev/stdout | \
      STAR \
      --genomeDir genome_dir \
      --outStd Log \
      --readFilesIn /dev/stdin \
      --runThreadN ~{cpu} \
      --outSAMtype BAM Unsorted \
      --outFileNamePrefix "~{sample_id}_" \
      ~{" " + flags}

        gsutil -q -m cp ~{sample_id}_* ~{output_directory}/
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil(size(genome_tar, "GB") * 5 + (disk_space_multiplier * size(input_bam, "GB")) + extra_disk_space) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
  }

}
task AddTags {
  input {
    String memory
    File aligned_bam
    File unaligned_bam
    String sample_id
    Int preemptible
    File genome_fasta
    File genome_dict
    String zones
    File refflat
    File gene_intervals
    String output_directory
    String drop_seq_tools_version
    Float disk_space_multiplier
    Int cpu
    Int sort_bam_max_records_in_ram
    String docker_registry
  }


  output {
    String bam = "${output_directory}/${sample_id}_aligned_tagged.bam"
    String output_sample_id = "${sample_id}"
  }
  command <<<

      set -o pipefail
       set -e

      mkfifo pipe1 pipe2 pipe3

      java -Dsamjdk.compression_level=1 -Xms3000m -jar /software/picard.jar SortSam VALIDATION_STRINGENCY=SILENT \
      INPUT=~{aligned_bam} \
      OUTPUT=pipe1 \
      MAX_RECORDS_IN_RAM=~{sort_bam_max_records_in_ram} \
      SORT_ORDER=queryname | \
      java -Dsamjdk.compression_level=1 -Xms3000m -jar /software/picard.jar MergeBamAlignment VALIDATION_STRINGENCY=SILENT \
      ALIGNED_BAM=pipe1 \
      UNMAPPED_BAM=~{unaligned_bam} \
      OUTPUT=pipe2 \
      REFERENCE_SEQUENCE=~{genome_fasta} \
      INCLUDE_SECONDARY_ALIGNMENTS=false \
      PAIRED_RUN=false \
      CLIP_ADAPTERS=false | \
      java -Dsamjdk.compression_level=1 -Xms3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagReadWithInterval VALIDATION_STRINGENCY=SILENT \
      I=pipe2 \
      O=pipe3 \
      INTERVALS=~{gene_intervals} \
      TAG=XG | \
      java -Dsamjdk.compression_level=2 -Xms3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagReadWithGeneFunction VALIDATION_STRINGENCY=SILENT \
      INPUT=pipe3 \
      O=~{sample_id}_aligned_tagged.bam \
      ANNOTATIONS_FILE=~{refflat}

      gsutil -q -m cp ~{sample_id}_aligned_tagged.bam ~{output_directory}/
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil(20 + size(genome_fasta, "GB") + size(gene_intervals, "GB") + size(genome_dict, "GB") + size(refflat, "GB") + size(unaligned_bam, "GB") * 2 + size(aligned_bam, "GB") * disk_space_multiplier) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
  }

}

