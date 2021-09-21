version 1.0


workflow dropest {
  input {
    Int max_umi_merge_edit_distance = 1
    String gene_id_tag = "XG"
    Boolean velocyto = true
    String sample_id
    File input_bam
    String zones = "us-east1-d us-west1-a us-west1-b"
    Boolean merge_barcodes_precise = true
    Float min_merge_fraction = 0.2
    String docker_registry
    String dropest_version = "0.8.5"
    String output_directory
    String cellular_barcode_quality_tag = "CY"
    String umi_quality_tag = "UY"
    String umi_tag = "XM"
    String dropest_memory = "104G"
    Int preemptible = 2
    String cellular_barcode_tag = "XC"
    Boolean apply_directional_umi_correction = true
    Int max_cb_merge_edit_distance = 2
    Int dropest_cpu = 1
    Int? genes_min = 100
    File? cellular_barcode_whitelist
    Int min_genes_before_merge = 20
    Int? cells_max
  }
  call run_dropest {
    input:
      memory = dropest_memory,
      input_bam = input_bam,
      sample_id = sample_id,
      apply_directional_umi_correction = apply_directional_umi_correction,
      merge_barcodes_precise = merge_barcodes_precise,
      cellular_barcode_whitelist = cellular_barcode_whitelist,
      genes_min = genes_min,
      cells_max = cells_max,
      min_merge_fraction = min_merge_fraction,
      max_cb_merge_edit_distance = max_cb_merge_edit_distance,
      max_umi_merge_edit_distance = max_umi_merge_edit_distance,
      min_genes_before_merge = min_genes_before_merge,
      cellular_barcode_tag = cellular_barcode_tag,
      umi_tag = umi_tag,
      gene_id_tag = gene_id_tag,
      cellular_barcode_quality_tag = cellular_barcode_quality_tag,
      umi_quality_tag = umi_quality_tag,
      velocyto = velocyto,
      cpu = dropest_cpu,
      preemptible = preemptible,
      zones = zones,
      output_directory = output_directory,
      dropest_version = dropest_version,
      docker_registry = docker_registry
  }

  output {
    String? count_matrices = run_dropest.count_matrices
    String count_matrix = run_dropest.count_matrix
    String bam = run_dropest.bam
    String log = run_dropest.log
  }
}
task run_dropest {
  input {
    String memory
    File input_bam
    String input_bam_name = basename(input_bam, ".bam")
    String sample_id
    Boolean apply_directional_umi_correction
    Boolean merge_barcodes_precise
    File? cellular_barcode_whitelist
    Int? genes_min
    Int? cells_max
    Float min_merge_fraction
    Int max_cb_merge_edit_distance
    Int max_umi_merge_edit_distance
    Int min_genes_before_merge
    String cellular_barcode_tag
    String umi_tag
    String gene_id_tag
    String cellular_barcode_quality_tag
    String umi_quality_tag
    Boolean velocyto
    Int cpu
    Int preemptible
    String zones
    String output_directory
    String dropest_version
    String docker_registry
  }


  output {
    String? count_matrices = "${output_directory}/${sample_id}.cell.counts.matrices.rds"
    String count_matrix = "${output_directory}/${sample_id}.cell.counts.rds"
    String bam = "${output_directory}/${input_bam_name}.dropest_filtered.bam"
    String log = "${output_directory}/est_main.log"
  }
  command <<<

      set -e

      create_dropest_config.py \
      --output dropest_config.xml \
      --min_merge_fraction ~{min_merge_fraction} \
        --max_cb_merge_edit_distance ~{max_cb_merge_edit_distance} \
        --max_umi_merge_edit_distance ~{max_umi_merge_edit_distance} \
        --min_genes_before_merge ~{min_genes_before_merge} \
        --cb ~{cellular_barcode_tag} \
        --umi ~{umi_tag} \
        --gene ~{gene_id_tag} \
        --cb_quality ~{cellular_barcode_quality_tag} \
        --umi_quality ~{umi_quality_tag} \
        ~{"--whitelist " + cellular_barcode_whitelist}

      dropest -F \
      ~{true="-V " false=""  velocyto} \
      ~{true="-u " false=""  apply_directional_umi_correction} \
      ~{true="-M " false=""  merge_barcodes_precise} \
      ~{"-C " + cells_max} \
      ~{"-G " + genes_min} \
      -c dropest_config.xml \
      -f ~{input_bam}

      mv ~{input_bam_name}.filtered.bam ~{input_bam_name}.dropest_filtered.bam
      mv cell.counts.rds ~{sample_id}.cell.counts.rds


      if [ '~{velocyto}' == 'true' ]; then
         mv cell.counts.matrices.rds ~{sample_id}.cell.counts.matrices.rds
      fi

        gsutil -q -m cp est_main.log *.rds *.bam ~{output_directory}/
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil((3 * size(input_bam, "GB")) + 20) + " HDD"
    docker: "${docker_registry}/dropest:${dropest_version}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
  }

}
task run_dropest_umi_correct {
  input {
    String memory
    File count_matrix
    Int cpu
    String sample_id
    Int preemptible
    String zones
    String output_directory
    String dropest_version
    String docker_registry
  }


  output {
    String umi_corrected_count_matrix = "${output_directory}/${sample_id}.umi_corrected.rds"
  }
  command <<<

      set -e

      monitor_script.sh &

      dropest_umi_correct.R ~{count_matrix} ~{cpu}

      mv umi_corrected.rds ~{sample_id}.umi_corrected.rds
        gsutil -q -m cp *.rds ~{output_directory}/
   
  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil((1.5 * size(count_matrix, "GB")) + 20) + " HDD"
    docker: "${docker_registry}/dropest:${dropest_version}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
  }

}

