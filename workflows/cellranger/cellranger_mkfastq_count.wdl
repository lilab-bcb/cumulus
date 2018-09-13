import "https://api.firecloud.org/ga4gh/v1/tools/regev:cellranger_mkfastq/versions/3/plain-WDL/descriptor" as crm
import "https://api.firecloud.org/ga4gh/v1/tools/regev:cellranger_count/versions/16/plain-WDL/descriptor" as crc
import "https://api.firecloud.org/ga4gh/v1/tools/regev:scrtools_adt/versions/2/plain-WDL/descriptor" as sa

workflow cellranger_mkfastq_count {
	# 5 or 6 columns (Sample, Reference, Flowcell, Lane, Index, [Chemistry]). gs URL
	File input_csv_file
	# CellRanger output directory, gs URL
	String cellranger_output_directory

	# If run cellranger mkfastq
	Boolean run_mkfastq = true
	# If run cellranger count
	Boolean run_count = true


	# Whether to delete input_bcl_directory. If false, you should delete this folder yourself so as to not incur storage charges.
	Boolean delete_input_bcl_directory = false

	# Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
	Boolean? secondary = false
	# If force cells, default: true
	Boolean? do_force_cells = false
	# Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells. Default: 6,000 cells
	Int? force_cells = 6000
	# Expected number of recovered cells. Default: 3,000 cells. Mutually exclusive with force_cells
	Int? expect_cells = 3000

	# For extracting ADT count

	# antibody barcodes in csv format
	File? antibody_barcode_file

	# maximum hamming distance in antibody barcodes
	Int? max_mismatch = 3



	# 2.1.1 or 2.2.0
	String? cellranger_version = "2.1.1"

	# Number of cpus per cellranger job
	Int? num_cpu = 64
	# Memory in GB
	Int? memory = 128
	# Optional memory in GB for scrtools_adt
	Int? adt_memory = 32
	# Optional disk space for mkfastq.
	Int? mkfastq_disk_space = 1500
	# Optional disk space needed for cell ranger count.
	Int? count_disk_space = 500	
	# Optional disk space needed for scrtools_adt
	Int? adt_disk_space = 100
	# Number of preemptible tries 
	Int? preemptible = 2


	if (run_mkfastq) {
		call generate_bcl_csv {
			input:
				input_csv_file = input_csv_file,
				output_dir = cellranger_output_directory,
				cellranger_version = cellranger_version,
				preemptible = preemptible
		}

		scatter (run_id in generate_bcl_csv.run_ids) {
			call crm.cellranger_mkfastq as cellranger_mkfastq {
				input:
					input_bcl_directory = generate_bcl_csv.inpdirs[run_id],
					input_csv_file = generate_bcl_csv.bcls[run_id],
					output_directory = cellranger_output_directory,
					delete_input_bcl_directory = delete_input_bcl_directory,
					cellranger_version = cellranger_version,
					num_cpu = num_cpu,
					memory = memory,
					disk_space = mkfastq_disk_space,
					preemptible = preemptible
			}
		}
	}

	if (run_count) {
		call generate_count_config {
			input:
				input_csv_file = input_csv_file,
				run_ids = generate_bcl_csv.run_ids,
				fastq_dirs = cellranger_mkfastq.output_fastqs_flowcell_directory,
				cellranger_version = cellranger_version,
				preemptible = preemptible			
		}

		if (generate_count_config.sample_ids[0] != '') {
			scatter (sample_id in generate_count_config.sample_ids) {
				call crc.cellranger_count as cellranger_count {
					input:
						sample_id = sample_id,
						input_fastqs_directories = generate_count_config.sample2dir[sample_id],
						output_directory = cellranger_output_directory,
						genome = generate_count_config.sample2genome[sample_id],
						chemistry = generate_count_config.sample2chemistry[sample_id],
						secondary = secondary,
						do_force_cells = do_force_cells,
						force_cells = force_cells,
						expect_cells = expect_cells,
						cellranger_version = cellranger_version,
						num_cpu = num_cpu,
						memory = memory,
						disk_space = count_disk_space,
						preemptible = preemptible
				}
			}

			call collect_summaries {
				input:
					summaries = cellranger_count.output_metrics_summary,
					sample_ids = cellranger_count.output_count_directory,
					cellranger_version = cellranger_version,
					preemptible = preemptible
			}		
		}

		if (generate_count_config.sample_adt_ids[0] != '' && defined(antibody_barcode_file)) {
			scatter (sample_id in generate_count_config.sample_adt_ids) {
				call sa.scrtools_adt as scrtools_adt {
					input:
						sample_id = sample_id,
						input_fastqs_directories = generate_count_config.sample_adt2dir[sample_id],
						output_directory = cellranger_output_directory,
						antibody_barcode_file = antibody_barcode_file,
						max_mismatch = max_mismatch,
						memory = adt_memory,
						disk_space = adt_disk_space,
						preemptible = preemptible
				}
			}		
		}
	}
}

task generate_bcl_csv {
	File input_csv_file
	String output_dir
	String cellranger_version
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import os
		import pandas as pd 
		from subprocess import check_call
		df = pd.read_csv('${input_csv_file}', header = 0)
		with open('run_ids.txt', 'w') as fo1, open('inpdirs.txt', 'w') as fo2, open('bcls.txt', 'w') as fo3:
			for input_dir in df['Flowcell'].unique():
				run_id = os.path.basename(input_dir)
				bcl_df = df.loc[df['Flowcell'] == input_dir, ['Lane', 'Sample', 'Index']]
				bcl_df.to_csv(run_id + '_bcl.csv', index = False)
				call_args = ['gsutil', '-q', 'cp', run_id + '_bcl.csv', '${output_dir}/']
				# call_args = ['cp', run_id + '_bcl.csv', '${output_dir}/']
				print(' '.join(call_args))
				check_call(call_args)
				fo1.write(run_id + '\n')
				fo2.write(run_id + '\t' + input_dir + '\n')
				fo3.write(run_id + '\t${output_dir}/' + run_id + '_bcl.csv\n')
		CODE
	}

	output {
		Array[String] run_ids = read_lines('run_ids.txt')
		Map[String, String] inpdirs = read_map('inpdirs.txt')
		Map[String, String] bcls = read_map('bcls.txt')
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		preemptible: "${preemptible}"
	}
}

task generate_count_config {
	File input_csv_file
	Array[String]? run_ids
	Array[String]? fastq_dirs
	String cellranger_version
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import os
		import pandas as pd
		from subprocess import check_call
		df = pd.read_csv('${input_csv_file}', header = 0)
		run_ids = '${sep="," run_ids}'.split(',')
		fastq_dirs = '${sep="," fastq_dirs}'.split(',')
		rid2fdir = dict()
		for run_id, fastq_dir in zip(run_ids, fastq_dirs):
			if run_id is not '':
				rid2fdir[run_id] = fastq_dir
		with open('sample_ids.txt', 'w') as fo1, open('sample2dir.txt', 'w') as fo2, open('sample2genome.txt', 'w') as fo3, open('sample2chemistry.txt', 'w') as fo4, open('sample_adt_ids.txt', 'w') as fo5, open('sample_adt2dir.txt', 'w') as fo6:
			n_normal = 0
			n_non_normal = 0
			for sample_id in df['Sample'].unique():
				df_local = df.loc[df['Sample'] == sample_id]
				assert df_local['Index'].unique().size == 1
				is_adt = df_local['Index'].iat[0].find('-') < 0
				if not is_adt:
					n_normal += 1
					fo1.write(sample_id + '\n')
					dirs = df_local['Flowcell'].map(lambda x: x if len(rid2fdir) == 0 else rid2fdir[os.path.basename(x)]).values
					fo2.write(sample_id + '\t' + ','.join(dirs) + '\n')
					assert df_local['Reference'].unique().size == 1
					fo3.write(sample_id + '\t' + df_local['Reference'].iat[0] + '\n')
					chemistry = 'auto'
					if 'Chemistry' in df_local.columns:
						assert df_local['Chemistry'].unique().size == 1
						chemistry = df_local['Chemistry'].iat[0]
					fo4.write(sample_id + '\t' + chemistry + '\n')
				else:
					n_non_normal += 1
					fo5.write(sample_id + '\n')
					dirs = df_local['Flowcell'].map(lambda x: x if len(rid2fdir) == 0 else rid2fdir[os.path.basename(x)]).values
					fo6.write(sample_id + '\t' + ','.join(dirs) + '\n')
			if n_normal == 0:
				fo2.write('null\tnull\n')
				fo3.write('null\tnull\n')
				fo4.write('null\tnull\n')
			if n_non_normal == 0:
				fo5.write('null\tnull\n')
				fo6.write('null\tnull\n')
		CODE
	}

	output {
		Array[String] sample_ids = read_lines('sample_ids.txt')
		Map[String, String] sample2dir = read_map('sample2dir.txt')
		Map[String, String] sample2genome = read_map('sample2genome.txt')
		Map[String, String] sample2chemistry = read_map('sample2chemistry.txt')
		Array[String] sample_adt_ids = read_lines('sample_adt_ids.txt')
		Map[String, String] sample_adt2dir = read_map('sample_adt2dir.txt')
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		preemptible: "${preemptible}"
	}
}

task collect_summaries {
	Array[File] summaries
	Array[String] sample_ids
	String? cellranger_version
	Int? preemptible

	command {
		python <<CODE
		import pandas as pd
		import os
		import xlsxwriter
		summaries = pd.read_csv('${write_lines(summaries)}', header = None)
		sample_ids = pd.read_csv('${write_lines(sample_ids)}', header = None).applymap(lambda x: os.path.basename(x))
		df_info = pd.concat([summaries, sample_ids], axis = 1)
		df_info.columns = ['summary', 'sample_id']
		dfs = []
		for idx, row in df_info.iterrows():
			df = pd.read_csv(row['summary'], header = 0)
			df.index = [row['sample_id']]
			dfs.append(df)
		df_res = pd.concat(dfs)
		df_res.index.name = "Sample"
		writer = pd.ExcelWriter('summaries.xlsx', engine = 'xlsxwriter')
		df_res.to_excel(writer, sheet_name = "summaries")
		writer.save()
		CODE
	}

	output {
		File metrics_summaries = "summaries.xlsx"
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		preemptible: "${preemptible}"
	}
}
