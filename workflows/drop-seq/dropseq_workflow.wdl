import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:bcl2fastq/versions/2/plain-WDL/descriptor" as bcl2fastq_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropest/versions/3/plain-WDL/descriptor" as dropest_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropseq_align/versions/4/plain-WDL/descriptor" as dropseq_align_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropseq_count/versions/4/plain-WDL/descriptor" as dropseq_count_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropseq_prepare_fastq/versions/4/plain-WDL/descriptor" as dropseq_prepare_fastq_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:dropseq_qc/versions/3/plain-WDL/descriptor" as dropseq_qc_wdl

workflow dropseq_workflow {
    # Either a list of flowcell URLs or a tab separated file with no header with sample_id tab r1 tab r2
    File input_tsv_file
    # Output directory, gs URL
    String output_directory
    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "/+$", "")

    Array[Array[String]] input_tsv = read_tsv(input_tsv_file)
    Boolean run_bcl2fastq = false
    Boolean run_dropseq_tools = true
    Boolean run_dropest = false

    String? star_flags

    # hg19, mm10, hg19_mm10, mmul_8.0.1 or a path to a custom reference JSON file
    String reference
    File? acronym_file = "gs://regev-lab/resources/DropSeq/index.json"
    String? bcl2fastq_docker_registry = "gcr.io/broad-cumulus"
    String bcl2fastq_docker_registry_stripped = sub(bcl2fastq_docker_registry, "/+$", "")
    String docker_registry = "cumulusprod"
    # docker_registry with trailing slashes stripped
    String docker_registry_stripped = sub(docker_registry, "/+$", "")

    # use ncells value directly instead of estimating from elbow plot
    Int? drop_seq_tools_force_cells
    File? cellular_barcode_whitelist

    # Minimal number of genes for cells after the merge procedure
    Int? dropest_genes_min = 100
    # maximal number of output cells
    Int? dropest_cells_max
    #  Threshold for the merge procedure
    Float dropest_min_merge_fraction  = 0.2
    # Max edit distance between barcodes.
    Int? dropest_max_cb_merge_edit_distance = 2

    # Max edit distance between UMIs
    Int? dropest_max_umi_merge_edit_distance = 1

    # Minimal number of genes for cells before the merge procedure. Used mostly for optimization.
    Int? dropest_min_genes_before_merge = 10

    # use precise merge strategy (can be slow), recommended to use when the list of real barcodes is not available
    Boolean? dropest_merge_barcodes_precise = true

    # save separate count matrices for exons, introns and exon/intron spanning reads
    Boolean? dropest_velocyto = true
    # Whether to delete input_bcl_directory
    Boolean? delete_input_bcl_directory = false
    # Number of preemptible tries
    Int? preemptible = 2

    String? umi_base_range="13-20"
    String? cellular_barcode_base_range = "1-12"

    # The sequence to look for at the start of reads for trimming
    String? trim_sequence = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"

    # How many bases at the begining of the sequence must match before trimming occurs
    Int? trim_num_bases = 5

    Float? prepare_fastq_disk_space_multiplier = 4
    Float? add_bam_tags_disk_space_multiplier = 25

    Int? star_cpus
    # specify memory to override default for genome
    String? star_memory

    # extra disk space for STAR
    Float star_extra_disk_space = 3
    # multiply size of input bam by this factor
    Float star_disk_space_multiplier = 4

    Int? bcl2fastq_cpu = 64
    String? bcl2fastq_memory = "57.6G"
    Int? bcl2fastq_disk_space = 1500
    Int bcl2fastq_minimum_trimmed_read_length = 10
    Int bcl2fastq_mask_short_adapter_reads = 10
    String? dropest_memory = "104G"
    String? zones = "us-east1-d us-west1-a us-west1-b"
    String? drop_seq_tools_version = "2.3.0"
    String? bcl2fastq_version = "2.20.0.422"
    String? dropseq_report_version = "1.0.0"
    String? dropest_version = "0.8.6"
    String? merge_bam_alignment_memory="13G"
    Int? sort_bam_max_records_in_ram = 2000000
    String? drop_deq_tools_prep_bam_memory = "3750M"
    String? drop_deq_tools_dge_memory = "3750M"
    Array[String]? species # for mixed species

    if (run_bcl2fastq) {
        scatter (row in input_tsv) {
            call bcl2fastq_wdl.bcl2fastq {
                input:
                    input_bcl_directory = row[0],
                    output_directory = output_directory_stripped,
                    delete_input_bcl_directory = delete_input_bcl_directory,
                    minimum_trimmed_read_length = bcl2fastq_minimum_trimmed_read_length,
                    mask_short_adapter_reads = bcl2fastq_mask_short_adapter_reads,
                    zones = zones,
                    num_cpu = bcl2fastq_cpu,
                    memory = bcl2fastq_memory,
                    disk_space = bcl2fastq_disk_space,
                    preemptible = preemptible,
                    bcl2fastq_version=bcl2fastq_version,
                    docker_registry=bcl2fastq_docker_registry_stripped
            }
        }

    }


    call generate_count_config {
        input:
            input_csv_file = input_tsv_file,
            drop_seq_tools_version=drop_seq_tools_version,
            bcl2fastq_sample_sheets = bcl2fastq.fastqs,
            disk_space_multiplier=prepare_fastq_disk_space_multiplier,
            zones = zones,
            preemptible = preemptible,
            acronym_file=acronym_file,
            star_cpus = star_cpus,
            star_memory = star_memory,
            species = species,
            reference=reference,
            docker_registry=docker_registry_stripped,
    }

    scatter (row in generate_count_config.grouped_sample_sheet) {
        call dropseq_prepare_fastq_wdl.dropseq_prepare_fastq as dropseq_prepare_fastq {
            input:
                r1 = row[1],
                r2 = row[2],
                disk_space = row[3],
                sample_id = row[0],
                umi_base_range=umi_base_range,
                trim_sequence = trim_sequence,
                trim_num_bases = trim_num_bases,
                cellular_barcode_base_range=cellular_barcode_base_range,
                quality_tags = drop_seq_tools_version != "2.1.0",
                drop_seq_tools_version=drop_seq_tools_version,
                output_directory = output_directory_stripped + '/' + row[0],
                zones = zones,
                preemptible = preemptible,
                docker_registry = docker_registry_stripped
        }


        call dropseq_align_wdl.dropseq_align as dropseq_align {
            input:
                sample_id = row[0],
                drop_seq_tools_version=drop_seq_tools_version,
                add_bam_tags_disk_space_multiplier=add_bam_tags_disk_space_multiplier,
                output_directory = output_directory_stripped + '/' + row[0],
                input_bam = dropseq_prepare_fastq.trimmed_bam,
                star_cpus = generate_count_config.star_cpus_output,
                star_memory = generate_count_config.star_memory_output,
                star_extra_disk_space = star_extra_disk_space,
                star_disk_space_multiplier = star_disk_space_multiplier,
                star_flags = star_flags,
                star_genome_file= generate_count_config.star_genome,
                refflat=generate_count_config.refflat,
                gene_intervals=generate_count_config.gene_intervals,
                genome_fasta=generate_count_config.genome_fasta,
                genome_dict=generate_count_config.genome_dict,
                merge_bam_alignment_memory=merge_bam_alignment_memory,
                sort_bam_max_records_in_ram =sort_bam_max_records_in_ram,
                zones = zones,
                preemptible = preemptible,
                docker_registry = docker_registry_stripped
        }

        if(run_dropest) {
            call dropest_wdl.dropest as dropest {
                input:
                    sample_id = row[0],
                    output_directory = output_directory_stripped + '/' + row[0],
                    input_bam = dropseq_align.aligned_tagged_bam,
                    velocyto=dropest_velocyto,
                    genes_min = dropest_genes_min,
                    cells_max = dropest_cells_max,
                    min_merge_fraction=dropest_min_merge_fraction,
                    max_cb_merge_edit_distance=dropest_max_cb_merge_edit_distance,
                    max_umi_merge_edit_distance=dropest_max_umi_merge_edit_distance,
                    min_genes_before_merge=dropest_min_genes_before_merge,
                    dropest_memory = dropest_memory,
                    merge_barcodes_precise = dropest_merge_barcodes_precise,
                    cellular_barcode_whitelist=cellular_barcode_whitelist,
                    dropest_version=dropest_version,
                    zones = zones,
                    preemptible = preemptible,
                    docker_registry=docker_registry_stripped
            }
        }

        if(run_dropseq_tools) {
            call dropseq_count_wdl.dropseq_count as dropseq_count {
                input:
                    sample_id = row[0],
                    dge_prep_memory = drop_deq_tools_prep_bam_memory,
                    dge_memory = drop_deq_tools_dge_memory,
                    drop_seq_tools_version=drop_seq_tools_version,
                    species=generate_count_config.species_list,
                    output_directory = output_directory_stripped + '/' + row[0],
                    input_bam = dropseq_align.aligned_tagged_bam,
                    force_cells = drop_seq_tools_force_cells,
                    cellular_barcode_whitelist=cellular_barcode_whitelist,
                    zones = zones,
                    preemptible = preemptible,
                    docker_registry=docker_registry_stripped
            }
        }


        call dropseq_qc_wdl.dropseq_qc as dropseq_qc {
            input:
                sample_id = row[0],
                input_bam = dropseq_align.aligned_tagged_bam,
                cell_barcodes=dropseq_count.cell_barcodes,
                refflat=generate_count_config.refflat,
                drop_seq_tools_version=drop_seq_tools_version,
                output_directory = output_directory_stripped + '/' + row[0],
                zones = zones,
                preemptible = preemptible,
                docker_registry=docker_registry_stripped
        }

    }


    call collect_summary {
        input:
            dge_summary=dropseq_count.dge_summary,
            dge_summary_multi_species=dropseq_count.dge_summary_multi_species,
            bead_synthesis_summary=dropseq_count.bead_synthesis_summary,
            star_log_final=dropseq_align.star_log_final,
            adapter_trimming_report=dropseq_prepare_fastq.adapter_trimming_report,
            polyA_trimming_report= dropseq_prepare_fastq.polyA_trimming_report,
            sc_rnaseq_metrics_report=dropseq_qc.sc_rnaseq_metrics_report,
            zones = zones,
            preemptible = preemptible,
            sample_id=dropseq_align.output_sample_id,
            output_directory=output_directory,
            docker_registry=docker_registry_stripped,
            version=dropseq_report_version
    }

    output {
        Array[String] trimmed_bam=dropseq_prepare_fastq.trimmed_bam
        Array[String] cellular_tag_summary=dropseq_prepare_fastq.cellular_tag_summary
        Array[String] molecular_tag_summary=dropseq_prepare_fastq.molecular_tag_summary
        Array[String] adapter_trimming_report = dropseq_prepare_fastq.adapter_trimming_report
        Array[String] polyA_trimming_report= dropseq_prepare_fastq.polyA_trimming_report

        Array[String] aligned_tagged_bam=dropseq_align.aligned_tagged_bam
        Array[String] aligned_bam=dropseq_align.aligned_bam
        Array[String] star_log_final = dropseq_align.star_log_final

        Array[String?] bead_synthesis_stats = dropseq_count.bead_synthesis_stats
        Array[String?] bead_synthesis_summary = dropseq_count.bead_synthesis_summary
        Array[String?] bead_synthesis_report = dropseq_count.bead_synthesis_report
        Array[String?] bead_substitution_report =dropseq_count.bead_substitution_report
        Array[String?] histogram=dropseq_count.histogram
        Array[String?] ncells=dropseq_count.ncells
        Array[String?] cumplot=dropseq_count.cumplot
        Array[String?] reads_plot=dropseq_count.reads_plot
        Array[String?] cell_barcodes = dropseq_count.cell_barcodes

        Array[String?] dge=dropseq_count.dge
        Array[String?] dge_summary = dropseq_count.dge_summary
        Array[String?] dge_reads=dropseq_count.dge_reads
        Array[String?] dge_summary_reads = dropseq_count.dge_summary_reads

        Array[Array[String]?] dge_multi_species = dropseq_count.dge_multi_species
        Array[Array[String]?] dge_summary_multi_species = dropseq_count.dge_summary_multi_species
        Array[Array[String]?] dge_reads_multi_species =dropseq_count.dge_reads_multi_species
        Array[Array[String]?] dge_summary_reads_multi_species = dropseq_count.dge_summary_reads_multi_species

        Array[String?] dropest_count_matrices = dropest.count_matrices
        Array[String?] dropest_count_matrix = dropest.count_matrix
        Array[String?] dropest_bam = dropest.bam
        Array[String?] dropest_log= dropest.log

        String sc_rnaseq_metrics_report=collect_summary.report
    }
}


task collect_summary {
    Array[String] sample_id
    Array[String?] dge_summary
    Array[Array[String]?] dge_summary_multi_species
    Array[String?] bead_synthesis_summary
    Array[String] star_log_final
    Array[String] adapter_trimming_report
    Array[String] polyA_trimming_report
    Array[String] sc_rnaseq_metrics_report

    String output_directory
    String zones
    Int preemptible
    String version
    String docker_registry

    command {
        set -e

        python /software/report.py \
        --sample_id '${sep="," sample_id}' \
        --dge_summary '${sep="," dge_summary}' \
        --dge_summary_multi_species '${sep="," dge_summary_multi_species}' \
        --star_log '${sep="," star_log_final}' \
        --adapter_trimming_report '${sep="," adapter_trimming_report}' \
        --polyA_trimming_report '${sep="," polyA_trimming_report}' \
        --sc_rnaseq_metrics_report '${sep="," sc_rnaseq_metrics_report}' \
        --bead_synthesis_summary '${sep="," bead_synthesis_summary}'

        gsutil -q -m cp drop_seq_report.html ${output_directory}/
    }

    output {
        String report = "${output_directory}/drop_seq_report.html"
    }

    runtime {
        cpu:1
        bootDiskSizeGb: 12
        disks: "local-disk 2 HDD"
        memory:"1GB"
        docker: "${docker_registry}/dropseq_report:${version}"
        zones: zones
        preemptible: "${preemptible}"
    }
}



task generate_count_config {
    File input_csv_file
    # array of sample_id tab r1,r1_n r2,r2_n
    Array[File]? bcl2fastq_sample_sheets
    String zones
    Int preemptible
    File acronym_file
    String drop_seq_tools_version
    String reference
    Int? star_cpus
    Array[String]? species
    String? star_memory
    Boolean is_reference_url = sub(reference, "^gs://.+", "URL") == "URL"
    File config_file = (if is_reference_url then reference else acronym_file)
    String docker_registry
    Float disk_space_multiplier

    command {
        set -e

        python <<CODE

        import json
        from subprocess import check_call

        import numpy as np
        import pandas as pd
        import os

        disk_space_multiplier = ${disk_space_multiplier}
        bcl2fastq_sample_sheets = '${sep="," bcl2fastq_sample_sheets}'.split(',')
        species = '${sep="," species}'.split(',')
        if len(species) == 1 and species[0] == '':
            species = []
        bcl2fastq_sample_sheets = list(filter(lambda x: x.strip() != '', bcl2fastq_sample_sheets))

        if len(bcl2fastq_sample_sheets) == 0:  # no bcl2fastq run, already a sample sheet with name, r1, r2
            df = pd.read_csv('${input_csv_file}', header=None, dtype=str, sep=None, engine='python')
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

        reference = '${reference}'.strip()
        with open('${config_file}', 'r') as f:
            config = json.load(f)
        if '${is_reference_url}' == 'false':
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
            star_memory = '${star_memory}'.strip()
            if star_memory != '':
                memory = star_memory
            w6.write(str(memory))
            cpu = config.get('star_cpus', '64')
            star_cpus = '${star_cpus}'.strip()
            if star_cpus != '':
                cpu = star_cpus
            w7.write(str(cpu))
            if len(species) == 0:
                species = config.get('species')
            if species is not None and len(species) > 0:
                for s in species:
                    w8.write(s + '\n')

        CODE
    }

    output {
        String star_genome = read_string('star_genome.txt')
        Array[String] species_list = read_lines('species.txt')
        String refflat = read_string('refflat.txt')
        String gene_intervals = read_string('gene_intervals.txt')
        String genome_fasta = read_string('genome_fasta.txt')
        String genome_dict = read_string('genome_dict.txt')
        String star_memory_output = read_string('star_memory.txt')
        Int star_cpus_output = read_int('star_cpu.txt')
        Array[Array[String]] grouped_sample_sheet = read_tsv('grouped_sample_sheet.txt')
    }

    runtime {
        cpu:1
        bootDiskSizeGb: 12
        disks: "local-disk 1 HDD"
        memory:"1GB"
        docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
        zones: zones
        preemptible: "${preemptible}"
    }
}
