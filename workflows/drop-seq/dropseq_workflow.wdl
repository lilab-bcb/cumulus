import "https://api.firecloud.org/ga4gh/v1/tools/regev:dropseq_align/versions/1/plain-WDL/descriptor" as dropseq_align_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/regev:dropseq_bcl2fastq/versions/1/plain-WDL/descriptor" as bcl2fastq_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/regev:dropseq_count/versions/1/plain-WDL/descriptor" as dropseq_count_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/regev:dropseq_prepare_fastq/versions/1/plain-WDL/descriptor" as dropseq_prepare_fastq_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/regev:dropseq_qc/versions/1/plain-WDL/descriptor" as dropseq_qc_wdl

workflow dropseq_workflow {
	# Either a list of flowcell URLS or sample_id tab r1 tab r2
	File input_csv_file
	# Output directory, gs URL
	String output_directory
	# Output directory, with trailing slashes stripped
	String output_directory_stripped = sub(output_directory, "/+$", "")
	Boolean run_bcl2fastq = true
	Array[Array[String]] input_tsv = read_tsv(input_csv_file)
	Boolean run_dropseq = true

	# for large genomes (e.g. human and mouse combined genome, use 96 cpus, "86.4G", or 32 cpus, 120G, or 64 cpus, 240GB - https://cloud.google.com/compute/docs/machine-types
	# STAR will take forever ($$$) without enough RAM
	Int star_cpus = 64
	String star_memory = "57.6G"
	String star_flags = "--limitOutSJcollapsed 1000000 --twopassMode Basic"
	File star_genome_file
	File refflat
	File gene_intervals
	File genome_fasta
	File genome_dict
	Int? bcl2fastq_cpu = 64
	String? bcl2fastq_memory = "57.6G"

	Int? bcl2fastq_disk_space = 1500
	String? zones = "us-east1-d us-west1-a us-west1-b"
	# use ncells value directly instead of estimating from elbow plot
	Int? force_cells

	String? workflow_version = "0.0.1"

	# Whether to delete input_bcl_directory, default: false
	Boolean? delete_input_bcl_directory = false
	# Number of preemptible tries
	Int? preemptible = 2
	String? bcl2fastq_version = "2.20.0.422-2"
	Int? add_bam_tags_disk_space_multiplier = 25

	if (run_bcl2fastq) {
		scatter (row in input_tsv) {
			call bcl2fastq_wdl.dropseq_bcl2fastq as bcl2fastq {
				input:
					input_bcl_directory = row[0],
					output_directory = output_directory_stripped,
					delete_input_bcl_directory = delete_input_bcl_directory,
					zones = zones,
					num_cpu = bcl2fastq_cpu,
					memory = bcl2fastq_memory,
					disk_space = bcl2fastq_disk_space,
					preemptible = preemptible,
					bcl2fastq_version=bcl2fastq_version
			}
		}

	}

	if(run_dropseq) {

		call generate_count_config {
			input:
				input_csv_file = input_csv_file,
				workflow_version=workflow_version,
				bcl2fastq_sample_sheets = bcl2fastq.fastqs,
				zones = zones,
				preemptible = preemptible
		}

		scatter (row in generate_count_config.grouped_sample_sheet) {
			call dropseq_prepare_fastq_wdl.dropseq_prepare_fastq as dropseq_prepare_fastq {
				input:
					r1 = row[1],
					r2 = row[2],
					disk_space = row[3],
					sample_id = row[0],
					workflow_version=workflow_version,
					output_directory = output_directory_stripped + '/' + row[0],
					zones = zones,
					preemptible = preemptible
			}


			call dropseq_align_wdl.dropseq_align as dropseq_align {
				input:
					sample_id = row[0],
					workflow_version=workflow_version,
					add_bam_tags_disk_space_multiplier=add_bam_tags_disk_space_multiplier,
					output_directory = output_directory_stripped + '/' + row[0],
					input_bam = dropseq_prepare_fastq.bam,
					star_cpus = star_cpus,
					star_memory = star_memory,
					star_flags = star_flags,
					star_genome_file=star_genome_file,
					refflat=refflat,
					gene_intervals=gene_intervals,
					genome_fasta=genome_fasta,
					genome_dict=genome_dict,
					zones = zones,
					preemptible = preemptible
			}

			call dropseq_count_wdl.dropseq_count as dropseq_count {
				input:
					sample_id = row[0],
					workflow_version=workflow_version,
					output_directory = output_directory_stripped + '/' + row[0],
					input_bam = dropseq_align.aligned_tagged_bam,
					force_cells = force_cells,
					zones = zones,
					preemptible = preemptible
			}

			call dropseq_qc_wdl.dropseq_qc as dropseq_qc {
				input:
					sample_id = row[0],
					input_bam = dropseq_align.aligned_tagged_bam,
					cell_barcodes=dropseq_count.cell_barcodes,
					refflat=refflat,
					workflow_version=workflow_version,
					output_directory = output_directory_stripped + '/' + row[0],
					zones = zones,
					preemptible = preemptible
			}

		}


		call collect_summary {
			input:
				dge_summary=dropseq_count.dge_summary,
				star_log_final=dropseq_align.star_log_final,
				adapter_trimming_report=dropseq_prepare_fastq.adapter_trimming_report,
				polyA_trimming_report= dropseq_prepare_fastq.polyA_trimming_report,
				sc_rnaseq_metrics_report=dropseq_qc.sc_rnaseq_metrics_report,
				bead_synthesis_summary=dropseq_count.bead_synthesis_summary,
				workflow_version=workflow_version,
				zones = zones,
				preemptible = preemptible,
				sample_id=dropseq_align.output_sample_id,
				output_directory=output_directory
		}

	}
}




task collect_summary {
	Array[String] sample_id
	Array[String] dge_summary
    Array[String] star_log_final
    Array[String] adapter_trimming_report
    Array[String] polyA_trimming_report
    Array[String] sc_rnaseq_metrics_report
    Array[String] bead_synthesis_summary
	String output_directory
	String zones
	Int preemptible
	String workflow_version

	command {
		set -e

		python /software/summary.py \
		--sample_id '${sep="," sample_id}' \
		--dge_summary '${sep="," dge_summary}' \
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
		docker: "regevlab/dropseq-${workflow_version}"
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
	String workflow_version

	command {
		set -e

		python <<CODE

		import pandas as pd
		import numpy as np
		from subprocess import check_call

		bcl2fastq_sample_sheets = '${sep="," bcl2fastq_sample_sheets}'.split(',')
		bcl2fastq_sample_sheets = list(filter(lambda x:x.strip() !='', bcl2fastq_sample_sheets))
			
		if len(bcl2fastq_sample_sheets) == 0: # no bcl2fastq run, already a sample sheet with name, r1, r2
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
				output_df = pd.read_csv('tmp.txt', sep=' ', engine='python')
				total = output_df.iloc[output_df.shape[0] - 1][0]
				total = total
				sizes.append(total)
			df[3] = sizes
			df[3] = np.ceil(1 + 4*(df[3]/1e9)).astype(int)
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
			df[3] = np.ceil(1 + 4*(df[3]/1e9)).astype(int)
			df.to_csv('grouped_sample_sheet.txt', index=True, header=False, sep='\t')

		CODE


	}

	output {
		Array[Array[String]] grouped_sample_sheet = read_tsv('grouped_sample_sheet.txt')
	}

	runtime {
		cpu:1
		bootDiskSizeGb: 12
		disks: "local-disk 1 HDD"
		memory:"1GB"
		docker: "regevlab/dropseq-${workflow_version}"
		zones: zones
		preemptible: "${preemptible}"
	}
}
