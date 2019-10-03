import "https://api.firecloud.org/ga4gh/v1/tools/scCloud:smartseq2_per_plate/versions/13/plain-WDL/descriptor" as ss2pp
# import "smartseq2_per_plate.wdl" as ss2pp

workflow smartseq2 {
	# 4 columns (Cell, Plate, Read1, and Read2). gs URL
	File input_csv_file
	# Output directory, gs URL
	String output_directory
	# Output directory, with trailing slashes stripped
	String output_directory_stripped = sub(output_directory, "/+$", "")
	# Reference to align reads against, GRCm38, GRCh38, or mm10
	String reference

	# smartseq2 version, default to "1.0.0"
	String? smartseq2_version = "1.0.0"
	# Google cloud zones, default to "us-east1-d us-west1-a us-west1-b"
	String? zones = "us-east1-d us-west1-a us-west1-b"
	# Number of cpus per job
	Int? num_cpu = 4
	# Memory string
	String? memory = "3.60G"
	# factor to multiply size of R1 and R2 by for RSEM
	Float? disk_space_multiplier = 10
	# Number of preemptible tries 
	Int? preemptible = 2

    Int? generate_count_matrix_disk_space = 10

	call parse_input_csv {
		input:
			input_csv_file = input_csv_file,
			output_directory = output_directory_stripped,
			smartseq2_version = smartseq2_version,
			zones = zones,
			preemptible = preemptible
	}

	scatter (plate_name in parse_input_csv.plate_names) {
		call ss2pp.smartseq2_per_plate {
			input:
				sample_sheet = parse_input_csv.pn2ss[plate_name],
				plate_name = plate_name,
				output_directory = output_directory_stripped,
				reference = reference,
				smartseq2_version = smartseq2_version,
				zones = zones,
				num_cpu = num_cpu,
				memory = memory,
				disk_space_multiplier = disk_space_multiplier,
				generate_count_matrix_disk_space=generate_count_matrix_disk_space,
				preemptible = preemptible
		}
	}

	output {
        Array[Array[File]] rsem_gene = smartseq2_per_plate.rsem_gene
        Array[Array[File]] rsem_isoform = smartseq2_per_plate.rsem_isoform
        Array[Array[File]] rsem_time = smartseq2_per_plate.rsem_time
        Array[Array[File]] rsem_cnt =smartseq2_per_plate.rsem_cnt
        Array[Array[File]] rsem_model = smartseq2_per_plate.rsem_model
        Array[Array[File]] rsem_theta = smartseq2_per_plate.rsem_theta

        Array[String] output_count_matrix = smartseq2_per_plate.output_count_matrix
        Array[String] output_qc_report = smartseq2_per_plate.output_qc_report
    }
}

task parse_input_csv {
	File input_csv_file
	String output_directory
	String smartseq2_version
	String zones
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import pandas as pd 
		from subprocess import check_call

		df = pd.read_csv('${input_csv_file}', header = 0, dtype=str, index_col = 0)
		with open('plate_names.txt', 'w') as fo1, open('pn2ss.txt', 'w') as fo2:
			for plate_name in df['Plate'].unique():
				plate_df = df.loc[df['Plate'] == plate_name, ['Read1', 'Read2']]
				plate_df.to_csv(plate_name + '_sample_sheet.csv')
				call_args = ['gsutil', '-q', 'cp', plate_name + '_sample_sheet.csv', '${output_directory}/']
				# call_args = ['cp', plate_name + '_sample_sheet.csv', '${output_directory}/']
				print(' '.join(call_args))
				check_call(call_args)
				fo1.write(plate_name + '\n')
				fo2.write(plate_name + '\t${output_directory}/' + plate_name + '_sample_sheet.csv\n')
		CODE
	}

	output {
		Array[String] plate_names = read_lines('plate_names.txt')
		Map[String, String] pn2ss = read_map('pn2ss.txt')
	}

	runtime {
		docker: "cumulusprod/smartseq2:${smartseq2_version}"
		zones: zones
		preemptible: "${preemptible}"
	}
}
