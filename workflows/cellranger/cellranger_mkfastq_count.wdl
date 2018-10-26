import "https://api.firecloud.org/ga4gh/v1/tools/regev:cellranger_mkfastq/versions/10/plain-WDL/descriptor" as crm
import "https://api.firecloud.org/ga4gh/v1/tools/regev:cellranger_count/versions/18/plain-WDL/descriptor" as crc
import "https://api.firecloud.org/ga4gh/v1/tools/regev:cellranger_vdj/versions/1/plain-WDL/descriptor" as crv
import "https://api.firecloud.org/ga4gh/v1/tools/scCloud:scCloud_adt/versions/1/plain-WDL/descriptor" as sa

# import "../cellranger/cellranger_mkfastq.wdl" as crm
# import "../cellranger/cellranger_count.wdl" as crc
# import "../cellranger/cellranger_vdj.wdl" as crv
# import "../scCloud/scCloud_adt.wdl" as sa

workflow cellranger_mkfastq_count {
	# 5 or 6 columns (Sample, Reference, Flowcell, Lane, Index, [Chemistry]). gs URL
	File input_csv_file
	# Output directory, gs URL
	String output_directory
	# Output directory, with trailing slashes stripped
	String output_directory_stripped = sub(output_directory, "/+$", "")

	# If run cellranger mkfastq
	Boolean run_mkfastq = true
	# If run cellranger count
	Boolean run_count = true


	# Whether to delete input_bcl_directory, default: true
	Boolean? delete_input_bcl_directory = true

	# Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
	Boolean? secondary = false
	# If force cells, default: false
	Boolean? do_force_cells = false
	# Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells. Default: 6,000 cells
	Int? force_cells = 6000
	# Expected number of recovered cells. Default: 3,000 cells. Mutually exclusive with force_cells
	Int? expect_cells = 3000

	# For vdj

	# For vdj, force pipeline to use this number of cells, bypassing the cell detection algorithm.
	Int? vdj_force_cells
	# Do not align reads to reference V(D)J sequences before de novo assembly. Default: false
	Boolean? vdj_denovo = false
	# Force the web summary HTML and metrics summary CSV to only report on a particular chain type. The accepted values are: auto for autodetection based on TR vs IG representation, TR for T cell receptors, IG for B cell receptors, all for all chain types.
	String? vdj_chain 

	# For extracting ADT count

	# antibody barcodes in csv format
	File? antibody_barcode_file
	# maximum hamming distance in antibody barcodes
	Int? max_mismatch = 3

	# 2.2.0 or 2.1.1
	String? cellranger_version = "2.2.0"
	# Number of cpus per cellranger job
	Int? num_cpu = 64
	# Memory in GB
	Int? memory = 128
	# Optional memory in GB for scCloud_adt
	Int? adt_memory = 32
	# Optional disk space for mkfastq.
	Int? mkfastq_disk_space = 1500
	# Optional disk space needed for cell ranger count.
	Int? count_disk_space = 500	
	# Optional disk space needed for cell ranger vdj.
	Int? vdj_disk_space = 500	
	# Optional disk space needed for scCloud_adt
	Int? adt_disk_space = 100
	# Number of preemptible tries 
	Int? preemptible = 2

	if (run_mkfastq) {
		call generate_bcl_csv {
			input:
				input_csv_file = input_csv_file,
				output_dir = output_directory_stripped,
				cellranger_version = cellranger_version,
				preemptible = preemptible
		}

		scatter (run_id in generate_bcl_csv.run_ids) {
			call crm.cellranger_mkfastq as cellranger_mkfastq {
				input:
					input_bcl_directory = generate_bcl_csv.inpdirs[run_id],
					input_csv_file = generate_bcl_csv.bcls[run_id],
					output_directory = output_directory_stripped,
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
				output_dir = output_directory_stripped,
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
						output_directory = output_directory_stripped,
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

		if (generate_count_config.sample_vdj_ids[0] != '') {
			scatter (sample_id in generate_count_config.sample_vdj_ids) {
				call crv.cellranger_vdj as cellranger_vdj {
					input:
						sample_id = sample_id,
						input_fastqs_directories = generate_count_config.sample2dir[sample_id],
						output_directory = output_directory_stripped,
						genome = generate_count_config.sample2genome[sample_id],
						force_cells = vdj_force_cells,
						denovo = vdj_denovo,
						chain = vdj_chain,
						cellranger_version = cellranger_version,
						num_cpu = num_cpu,
						memory = memory,
						disk_space = vdj_disk_space,
						preemptible = preemptible
				}
			}

			call collect_summaries as collect_summaries_vdj {
				input:
					summaries = cellranger_vdj.output_metrics_summary,
					sample_ids = cellranger_vdj.output_vdj_directory,
					cellranger_version = cellranger_version,
					preemptible = preemptible
			}		
		}

		if (generate_count_config.sample_adt_ids[0] != '' && defined(antibody_barcode_file)) {
			scatter (sample_id in generate_count_config.sample_adt_ids) {
				call sa.scCloud_adt as scCloud_adt {
					input:
						sample_id = sample_id,
						input_fastqs_directories = generate_count_config.sample2dir[sample_id],
						output_directory = output_directory_stripped,
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
		import re
		import pandas as pd 
		from subprocess import check_call

		df = pd.read_csv('${input_csv_file}', header = 0)
		df['Flowcell'] = df['Flowcell'].map(lambda x: re.sub('/+$', '', x)) # remove trailing slashes
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
	String output_dir
	Array[String]? run_ids
	Array[String]? fastq_dirs
	String cellranger_version
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import os
		import re
		import pandas as pd
		from subprocess import check_call

		df = pd.read_csv('${input_csv_file}', header = 0)
		df['Sample'] = df['Sample'].astype(str)
		df['Flowcell'] = df['Flowcell'].map(lambda x: re.sub('/+$', '', x)) # remove trailing slashes

		run_ids = '${sep="," run_ids}'.split(',')
		fastq_dirs = '${sep="," fastq_dirs}'.split(',')
		rid2fdir = dict()
		for run_id, fastq_dir in zip(run_ids, fastq_dirs):
			if run_id is not '':
				rid2fdir[run_id] = fastq_dir

		with open('sample_ids.txt', 'w') as fo1, open('sample2dir.txt', 'w') as fo2, open('sample2genome.txt', 'w') as fo3, open('sample2chemistry.txt', 'w') as fo4, \
			 open('count_matrix.csv', 'w') as fo5, open('sample_vdj_ids.txt', 'w') as fo6, open('sample_adt_ids.txt', 'w') as fo7:

			fo5.write('Sample,Location\n')

			n_ref = n_chem = 0

			for sample_id in df['Sample'].unique():
				if sample_id.find(' ') != -1:
					raise ValueError('Invalid sample id: ' + sample_id)
				df_local = df.loc[df['Sample'] == sample_id]
				
				data_type = 'count'
				if 'DataType' in df_local.columns:
					assert df_local['DataType'].unique().size == 1
					data_type = df_local['DataType'].iat[0]

				if data_type == 'count':
					fo1.write(sample_id + '\n')
				elif data_type == 'vdj':
					fo6.write(sample_id + '\n')
				elif data_type == 'adt':
					fo7.write(sample_id + '\n')
				else:
					print('Invalid data type: ' + data_type + '!')
					assert False

				dirs = df_local['Flowcell'].map(lambda x: x if len(rid2fdir) == 0 else rid2fdir[os.path.basename(x)]).values
				fo2.write(sample_id + '\t' + ','.join(dirs) + '\n')
				
				if data_type != 'adt':
					assert df_local['Reference'].unique().size == 1
					fo3.write(sample_id + '\t' + df_local['Reference'].iat[0] + '\n')
					n_ref += 1

				if data_type == 'count':
					chemistry = 'auto'
					if 'Chemistry' in df_local.columns:
						assert df_local['Chemistry'].unique().size == 1
						chemistry = df_local['Chemistry'].iat[0]
					fo4.write(sample_id + '\t' + chemistry + '\n')
					n_chem += 1
					fo5.write(sample_id + ',${output_dir}/' + sample_id + '/filtered_gene_bc_matrices_h5.h5\n')

			if n_ref == 0:
				fo3.write('null\tnull\n')
			if n_chem == 0:
				fo4.write('null\tnull\n')
		CODE

		gsutil -q -m cp count_matrix.csv ${output_dir}/
		# cp count_matrix.csv ${output_dir}/
	}

	output {
		Array[String] sample_ids = read_lines('sample_ids.txt')
		Map[String, String] sample2dir = read_map('sample2dir.txt')
		Map[String, String] sample2genome = read_map('sample2genome.txt')
		Map[String, String] sample2chemistry = read_map('sample2chemistry.txt')
		String count_matrix = "${output_dir}/count_matrix.csv"
		Array[String] sample_vdj_ids = read_lines('sample_vdj_ids.txt')
		Array[String] sample_adt_ids = read_lines('sample_adt_ids.txt')
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
