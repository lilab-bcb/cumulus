import "https://api.firecloud.org/ga4gh/v1/tools/scCloud:smartseq2_per_plate/versions/2/plain-WDL/descriptor" as ss2pp
# import "../smartseq2/smartseq2_per_plate.wdl" as ss2pp

workflow smartseq2 {
	# 4 columns (Cell, Plate, Read1, and Read2). gs URL
	File input_csv_file
	# Output directory, gs URL
	String output_directory
	# Output directory, with trailing slashes stripped
	String output_directory_stripped = sub(output_directory, "/+$", "")
	# Reference to align reads against, now only GRCm38 is available
	String reference = "gs://regev-lab/resources/SmartSeq2/rsem_bowtie2_grcm38.tar.gz"
	# String reference = "rsem_bowtie2_grcm38.tar.gz"

	# Number of cpus per job
	Int? num_cpu = 4
	# Memory in GB
	Int? memory = 10
	# disk space in GB
	Int? disk_space = 10
	# Number of preemptible tries 
	Int? preemptible = 2

	call parse_input_csv {
		input:
			input_csv_file = input_csv_file,
			output_directory = output_directory_stripped,
			preemptible = preemptible
	}

	scatter (plate_name in parse_input_csv.plate_names) {
		call ss2pp.smartseq2_per_plate {
			input:
				sample_sheet = parse_input_csv.pn2ss[plate_name],
				plate_name = plate_name,
				output_directory = output_directory_stripped,
				reference = reference,
				num_cpu = num_cpu,
				memory = memory,
				disk_space = disk_space,
				preemptible = preemptible
		}
	}
}

task parse_input_csv {
	File input_csv_file
	String output_directory
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import pandas as pd 
		from subprocess import check_call

		df = pd.read_csv('${input_csv_file}', header = 0, index_col = 0)
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
		docker: "regevlab/smartseq2"
		preemptible: "${preemptible}"
	}
}
