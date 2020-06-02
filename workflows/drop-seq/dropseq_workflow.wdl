version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:bcl2fastq/versions/5/plain-WDL/descriptor" as bcl2fastq_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropest/versions/4/plain-WDL/descriptor"as dropest_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropseq_align/versions/5/plain-WDL/descriptor" as dropseq_align_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropseq_count/versions/5/plain-WDL/descriptor" as dropseq_count_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropseq_prepare_fastq/versions/5/plain-WDL/descriptor" as dropseq_prepare_fastq_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropseq_qc/versions/4/plain-WDL/descriptor" as dropseq_qc_wdl


# import "../bcl2fastq/bcl2fastq.wdl" as bcl2fastq_wdl
# import "dropest.wdl" as dropest_wdl
# import "dropseq_align.wdl" as dropseq_align_wdl
# import "dropseq_count.wdl" as dropseq_count_wdl
# import "dropseq_prepare_fastq.wdl" as dropseq_prepare_fastq_wdl
# import "dropseq_qc.wdl" as dropseq_qc_wdl

workflow dropseq_workflow {
  input {
    Float prepare_fastq_disk_space_multiplier = 4
    String drop_deq_tools_dge_memory = "3750M"
    Int sort_bam_max_records_in_ram = 2000000
    Int dropest_max_cb_merge_edit_distance = 2
    Boolean dropest_velocyto = true
    String merge_bam_alignment_memory = "13G"
    String dropest_version = "0.8.6"
    String reference
    String output_directory
    Int trim_num_bases = 5
    Int? dropest_cells_max
    String? star_flags
    File acronym_file = "gs://regev-lab/resources/DropSeq/index.json"
    Array[String]? species
    String bcl2fastq_version = "2.20.0.422"
    File? cellular_barcode_whitelist
    Boolean run_bcl2fastq = false
    String bcl2fastq_memory = "57.6G"
    String docker_registry = "cumulusprod"
    Int dropest_max_umi_merge_edit_distance = 1
    String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? dropest_genes_min = 100
    Boolean dropest_merge_barcodes_precise = true
    Float add_bam_tags_disk_space_multiplier = 25
    String? star_memory
    Int bcl2fastq_disk_space = 1500
    Float star_disk_space_multiplier = 4
    Int bcl2fastq_minimum_trimmed_read_length = 10
    Int preemptible = 2
    String dropest_memory = "104G"
    Int bcl2fastq_mask_short_adapter_reads = 10
    String umi_base_range = "13-20"
    File input_tsv_file
    Boolean run_dropest = false
    Boolean delete_input_bcl_directory = false
    Float star_extra_disk_space = 3
    String trim_sequence = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
    String drop_deq_tools_prep_bam_memory = "3750M"
    String dropseq_report_version = "1.0.1"
    Int? drop_seq_tools_force_cells
    Int bcl2fastq_cpu = 64
    String bcl2fastq_docker_registry = "gcr.io/broad-cumulus"
    String cellular_barcode_base_range = "1-12"
    Int? star_cpus
    Float dropest_min_merge_fraction = 0.2
    Int dropest_min_genes_before_merge = 10
    String drop_seq_tools_version = "2.3.0"
  }
  Array[Array[String]] input_tsv = read_tsv(input_tsv_file)
  String docker_registry_stripped = sub(docker_registry, "/+$", "")
  String bcl2fastq_docker_registry_stripped = sub(bcl2fastq_docker_registry, "/+$", "")
  String output_directory_stripped = sub(output_directory, "/+$", "")
  if (run_bcl2fastq) {
    scatter (row in input_tsv) {
      call bcl2fastq_wdl.bcl2fastq as bcl2fastq {
        input:
          preemptible = preemptible,
          delete_input_bcl_directory = delete_input_bcl_directory,
          memory = bcl2fastq_memory,
          output_directory = output_directory_stripped,
          num_cpu = bcl2fastq_cpu,
          zones = zones,
          minimum_trimmed_read_length = bcl2fastq_minimum_trimmed_read_length,
          bcl2fastq_version = bcl2fastq_version,
          input_bcl_directory = row[0],
          mask_short_adapter_reads = bcl2fastq_mask_short_adapter_reads,
          disk_space = bcl2fastq_disk_space,
          docker_registry = bcl2fastq_docker_registry_stripped
      }
    }
  }
  scatter (row in generate_count_config.grouped_sample_sheet) {
    if (run_dropest) {
      call dropest_wdl.dropest as dropest {
        input:
          min_merge_fraction = dropest_min_merge_fraction,
          zones = zones,
          velocyto = dropest_velocyto,
          max_umi_merge_edit_distance = dropest_max_umi_merge_edit_distance,
          min_genes_before_merge = dropest_min_genes_before_merge,
          docker_registry = docker_registry_stripped,
          preemptible = preemptible,
          output_directory = output_directory_stripped + "/" + row[0],
          cells_max = dropest_cells_max,
          merge_barcodes_precise = dropest_merge_barcodes_precise,
          cellular_barcode_whitelist = cellular_barcode_whitelist,
          dropest_memory = dropest_memory,
          input_bam = dropseq_align.aligned_tagged_bam,
          genes_min = dropest_genes_min,
          max_cb_merge_edit_distance = dropest_max_cb_merge_edit_distance,
          sample_id = row[0],
          dropest_version = dropest_version
      }
    }
    call dropseq_qc_wdl.dropseq_qc as dropseq_qc {
      input:
        drop_seq_tools_version = drop_seq_tools_version,
        docker_registry = docker_registry_stripped,
        output_directory = output_directory_stripped + "/" + row[0],
        zones = zones,
        preemptible = preemptible,
        refflat = generate_count_config.refflat,
        cell_barcodes = dropseq_count.cell_barcodes,
        input_bam = dropseq_align.aligned_tagged_bam,
        sample_id = row[0]
    }

  call dropseq_count_wdl.dropseq_count as dropseq_count {
    input:
      force_cells = drop_seq_tools_force_cells,
      drop_seq_tools_version = drop_seq_tools_version,
      docker_registry = docker_registry_stripped,
      dge_prep_memory = drop_deq_tools_prep_bam_memory,
      output_directory = output_directory_stripped + "/" + row[0],
      cellular_barcode_whitelist = cellular_barcode_whitelist,
      dge_memory = drop_deq_tools_dge_memory,
      zones = zones,
      input_bam = dropseq_align.aligned_tagged_bam,
      sample_id = row[0],
      species = generate_count_config.species_list,
      preemptible = preemptible
  }

    call dropseq_align_wdl.dropseq_align as dropseq_align {
      input:
        zones = zones,
        gene_intervals = generate_count_config.gene_intervals,
        drop_seq_tools_version = drop_seq_tools_version,
        docker_registry = docker_registry_stripped,
        output_directory = output_directory_stripped + "/" + row[0],
        preemptible = preemptible,
        star_flags = star_flags,
        sort_bam_max_records_in_ram = sort_bam_max_records_in_ram,
        star_extra_disk_space = star_extra_disk_space,
        merge_bam_alignment_memory = merge_bam_alignment_memory,
        star_genome_file = generate_count_config.star_genome,
        add_bam_tags_disk_space_multiplier = add_bam_tags_disk_space_multiplier,
        star_disk_space_multiplier = star_disk_space_multiplier,
        star_cpus = generate_count_config.star_cpus_output,
        refflat = generate_count_config.refflat,
        input_bam = dropseq_prepare_fastq.trimmed_bam,
        sample_id = row[0],
        genome_fasta = generate_count_config.genome_fasta,
        genome_dict = generate_count_config.genome_dict,
        star_memory = generate_count_config.star_memory_output
    }
    call dropseq_prepare_fastq_wdl.dropseq_prepare_fastq as dropseq_prepare_fastq {
      input:
        drop_seq_tools_version = drop_seq_tools_version,
        trim_num_bases = trim_num_bases,
        docker_registry = docker_registry_stripped,
        output_directory = output_directory_stripped + "/" + row[0],
        quality_tags = drop_seq_tools_version != "2.1.0",
        umi_base_range = umi_base_range,
        r2 = row[2],
        disk_space = row[3],
        preemptible = preemptible,
        cellular_barcode_base_range = cellular_barcode_base_range,
        r1 = row[1],
        trim_sequence = trim_sequence,
        sample_id = row[0],
        zones = zones
    }
  }

  call collect_summary {
    input:
      dge_summary = select_all(dropseq_count.dge_summary),
      dge_summary_multi_species = select_all(dropseq_count.dge_summary_multi_species),
      bead_synthesis_summary = dropseq_count.bead_synthesis_summary,
      sample_id = dropseq_align.output_sample_id,
      star_log_final = dropseq_align.star_log_final,
      adapter_trimming_report = dropseq_prepare_fastq.adapter_trimming_report,
      polyA_trimming_report = dropseq_prepare_fastq.polyA_trimming_report,
      sc_rnaseq_metrics_report = dropseq_qc.sc_rnaseq_metrics_report,
      output_directory = output_directory,
      zones = zones,
      preemptible = preemptible,
      version = dropseq_report_version,
      docker_registry = docker_registry_stripped
  }

  call generate_count_config {
    input:
      input_csv_file = input_tsv_file,
      bcl2fastq_sample_sheets = bcl2fastq.fastqs,
      zones = zones,
      preemptible = preemptible,
      acronym_file = acronym_file,
      drop_seq_tools_version = drop_seq_tools_version,
      reference = reference,
      star_cpus = star_cpus,
      species = species,
      star_memory = star_memory,
      docker_registry = docker_registry_stripped,
      disk_space_multiplier = prepare_fastq_disk_space_multiplier
  }

  output {
    Array[String?] dropest_count_matrices = dropest.count_matrices
    Array[String?] dge_summary = dropseq_count.dge_summary
    Array[String] aligned_bam = dropseq_align.aligned_bam
    Array[String] cellular_tag_summary = dropseq_prepare_fastq.cellular_tag_summary
    Array[Array[String]?] dge_summary_multi_species = dropseq_count.dge_summary_multi_species
    Array[String] star_log_final = dropseq_align.star_log_final
    Array[String?] bead_synthesis_summary = dropseq_count.bead_synthesis_summary
    Array[String?] bead_synthesis_stats = dropseq_count.bead_synthesis_stats
    Array[String] polyA_trimming_report = dropseq_prepare_fastq.polyA_trimming_report
    Array[String] molecular_tag_summary = dropseq_prepare_fastq.molecular_tag_summary
    Array[String?] cell_barcodes = dropseq_count.cell_barcodes
    Array[String?] ncells = dropseq_count.ncells
    Array[String?] reads_plot = dropseq_count.reads_plot
    Array[String?] bead_substitution_report = dropseq_count.bead_substitution_report
    Array[String?] dropest_bam = dropest.bam
    Array[String] aligned_tagged_bam = dropseq_align.aligned_tagged_bam
    String? sc_rnaseq_metrics_report = collect_summary.report
    Array[String?] dge_summary_reads = dropseq_count.dge_summary_reads
    Array[String?] histogram = dropseq_count.histogram
    Array[String?] dge = dropseq_count.dge
    Array[String] trimmed_bam = dropseq_prepare_fastq.trimmed_bam
    Array[String?] dropest_count_matrix = dropest.count_matrix
    Array[Array[String]?] dge_multi_species = dropseq_count.dge_multi_species
    Array[String?] bead_synthesis_report = dropseq_count.bead_synthesis_report
    Array[Array[String]?] dge_summary_reads_multi_species = dropseq_count.dge_summary_reads_multi_species
    Array[String] adapter_trimming_report = dropseq_prepare_fastq.adapter_trimming_report
    Array[String?] cumplot = dropseq_count.cumplot
    Array[String?] dge_reads = dropseq_count.dge_reads
    Array[Array[String]?] dge_reads_multi_species = dropseq_count.dge_reads_multi_species
    Array[String?] dropest_log = dropest.log
  }
}
task collect_summary {
  input {
    Array[String] dge_summary
    Array[Array[String]] dge_summary_multi_species
    Array[String] bead_synthesis_summary
    Array[String] sample_id
    Array[String] star_log_final
    Array[String] adapter_trimming_report
    Array[String] polyA_trimming_report
    Array[String] sc_rnaseq_metrics_report
    String output_directory
    String zones
    Int preemptible
    String version
    String docker_registry
  }


  output {
    String report = "${output_directory}/drop_seq_report.html"
  }
  command <<<

        set -e

        python /software/report.py \
        --sample_id '~{sep="," sample_id}' \
        --dge_summary '~{sep="," dge_summary}' \
        --dge_summary_multi_species ~{write_tsv(dge_summary_multi_species)} \
        --star_log '~{sep="," star_log_final}' \
        --adapter_trimming_report '~{sep="," adapter_trimming_report}' \
        --polyA_trimming_report '~{sep="," polyA_trimming_report}' \
        --sc_rnaseq_metrics_report '~{sep="," sc_rnaseq_metrics_report}' \
        --bead_synthesis_summary '~{sep="," bead_synthesis_summary}'

        gsutil -q -m cp drop_seq_report.html ~{output_directory}/

  >>>
  runtime {
    preemptible: "${preemptible}"
    bootDiskSizeGb: 12
    disks: "local-disk 2 HDD"
    docker: "${docker_registry}/dropseq_report:${version}"
    cpu: 1
    zones: zones
    memory: "1GB"
  }

}
task generate_count_config {
  input {
    File input_csv_file
    Array[File]? bcl2fastq_sample_sheets
    String zones
    Int preemptible
    File acronym_file
    String drop_seq_tools_version
    String reference
    Int? star_cpus
    Array[String]? species
    String? star_memory
    String docker_registry
    Float disk_space_multiplier
  }
    Boolean is_reference_url = sub(reference, "^gs://.+", "URL") == "URL"
    File config_file = (if is_reference_url then reference else acronym_file)


  output {
    String star_genome = read_string("star_genome.txt")
    Array[String] species_list = read_lines("species.txt")
    String refflat = read_string("refflat.txt")
    String gene_intervals = read_string("gene_intervals.txt")
    String genome_fasta = read_string("genome_fasta.txt")
    String genome_dict = read_string("genome_dict.txt")
    String star_memory_output = read_string("star_memory.txt")
    Int star_cpus_output = read_int("star_cpu.txt")
    Array[Array[String]] grouped_sample_sheet = read_tsv("grouped_sample_sheet.txt")
  }
  command <<<

        set -e

        python <<CODE

        import json
        from subprocess import check_call

        import numpy as np
        import pandas as pd
        import os

        disk_space_multiplier = ~{disk_space_multiplier}
        bcl2fastq_sample_sheets = '~{sep=","  bcl2fastq_sample_sheets}'.split(',')
        species = '~{sep=","  species}'.split(',')
        if len(species) == 1 and species[0] == '':
            species = []
        bcl2fastq_sample_sheets = list(filter(lambda x: x.strip() != '', bcl2fastq_sample_sheets))

        if len(bcl2fastq_sample_sheets) == 0:  # no bcl2fastq run, already a sample sheet with name, r1, r2
            df = pd.read_csv('~{input_csv_file}', header=None, dtype=str, sep=None, engine='python')
            for c in df.columns:
                df[c] = df[c].str.strip()
            df = df[[0, 1, 2]]
            df = df.groupby(0).agg(lambda col: ','.join(col))
            # get size for each sample fastqs
            sizes = []
            for i in range(df.shape[0]):
                disk_size_call = ['gsutil', 'du', '-c']
                disk_size_call = disk_size_call + df.iloc[i][1].split(',')
                disk_size_call = disk_size_call + df.iloc[i][2].split(',')
                output = open('tmp.txt', 'w')
                check_call(disk_size_call, stdout=output)
                output.close()
                output_df = pd.read_csv('tmp.txt', sep=' ', engine='python')
                if output_df.shape[0] == 0:
                    print(df.index[i] + ' fastqs not found')
                    exit(1)
                total = output_df.iloc[output_df.shape[0] - 1][0]
                total = total
                sizes.append(total)
            df[3] = sizes
            df[3] = np.ceil(1 + disk_space_multiplier * (df[3] / 1e9)).astype(int)
            df.to_csv('grouped_sample_sheet.txt', index=True, header=False, sep='\t')
        else:
            df = None
            for f in bcl2fastq_sample_sheets:
                df_i = pd.read_csv(f, header=None, dtype=str, sep='\t')
                df = pd.concat((df, df_i), copy=False) if df is not None else df_i
            agg_dict = dict()
            agg_dict[1] = lambda col: ','.join(col)
            agg_dict[2] = lambda col: ','.join(col)
            agg_dict[3] = 'sum'
            df = df.groupby(0).agg(agg_dict)
            df[3] = np.ceil(1 + disk_space_multiplier * (df[3].astype(float) / 1e9)).astype(int)
            df.to_csv('grouped_sample_sheet.txt', index=True, header=False, sep='\t')

        reference = '~{reference}'.strip()
        with open('~{config_file}', 'r') as f:
            config = json.load(f)
        if '~{is_reference_url}' == 'false':
            config = config[reference.lower()]
        with open('star_genome.txt', 'wt') as w1, open('refflat.txt', 'wt') as w2, open('gene_intervals.txt', 'wt') as w3, \
            open('genome_fasta.txt', 'wt') as w4, open('genome_dict.txt', 'wt') as w5, open('star_memory.txt', 'wt') as w6, \
            open('star_cpu.txt', 'wt') as w7, open('species.txt', 'wt') as w8:
            w1.write(config['star_genome'])
            w2.write(config['refflat'])
            w3.write(config['gene_intervals'])
            w4.write(config['genome_fasta'])
            w5.write(config['genome_dict'])
            memory = config.get('star_memory', '57.6G')
            star_memory = '~{star_memory}'.strip()
            if star_memory != '':
                memory = star_memory
            w6.write(str(memory))
            cpu = config.get('star_cpus', '64')
            star_cpus = '~{star_cpus}'.strip()
            if star_cpus != '':
                cpu = star_cpus
            w7.write(str(cpu))
            if len(species) == 0:
                species = config.get('species')
            if species is not None and len(species) > 0:
                for s in species:
                    w8.write(s + '\n')

        CODE

  >>>
  runtime {
    preemptible: "${preemptible}"
    bootDiskSizeGb: 12
    disks: "local-disk 1 HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "1GB"
  }

}
