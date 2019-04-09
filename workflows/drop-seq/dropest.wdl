workflow dropest {
   	String sample_id
    File input_bam
	String? zones = "us-east1-d us-west1-a us-west1-b"
	Int preemptible = 2
	String output_directory
	String dropest_version = "0.8.5"

	# minimal number of genes in output cells
	Int? genes_min
	# maximal number of output cells
	Int? cells_max
	String dropest_memory
	Int dropest_cpu = 1
#	String dropest_umi_correct_memory
#   	Int dropest_umi_correct_cpu
   	File? cellular_barcode_whitelist
	Boolean apply_directional_umi_correction = true
	# use precise merge strategy (can be slow), recommended to use when the list of real barcodes is not available
	Boolean merge_barcodes_precise = true
	call run_dropest {
		input:
			memory=dropest_memory,
			cpu=dropest_cpu,
			sample_id=sample_id,
			input_bam=input_bam,
			genes_min=genes_min,
			cells_max=cells_max,
			cellular_barcode_whitelist=cellular_barcode_whitelist,
			apply_directional_umi_correction=apply_directional_umi_correction,
			merge_barcodes_precise=merge_barcodes_precise,
			preemptible=preemptible,
			zones=zones,
			output_directory=output_directory,
			dropest_version=dropest_version
	}
#
#	if(!apply_directional_umi_correction){
#		call run_dropest_umi_correct {
#			input:
#				memory=dropest_umi_correct_memory,
#				cpu=dropest_umi_correct_cpu,
#				sample_id=sample_id,
#				count_matrix = run_dropest.count_matrix,
#				preemptible=preemptible,
#				zones=zones,
#				output_directory=output_directory,
#				dropest_version=dropest_version
#			}
#	}


	output {
		String count_matrices = run_dropest.count_matrices
		String count_matrix = run_dropest.count_matrix
		String log= run_dropest.log
#		String? umi_corrected_count_matrix = run_dropest_umi_correct.umi_corrected_count_matrix

	}

}



task run_dropest {
	String memory
	File input_bam
	Int cpu
	String sample_id
	Int preemptible
	String zones
	String output_directory
	String dropest_version
	Int? genes_min
	Int? cells_max
	File? cellular_barcode_whitelist
	Boolean apply_directional_umi_correction
    Boolean merge_barcodes_precise
	command {
		set -e

		create_dropest_config.py ${"--whitelist " + cellular_barcode_whitelist} --output dropest_config.xml
		dropest -V ${true="-u " false='' apply_directional_umi_correction} ${true="-M " false='' merge_barcodes_precise} ${"-C " + cells_max} ${"-G " + genes_min} -c dropest_config.xml -f ${input_bam}

		mv cell.counts.rds ${sample_id}.cell.counts.rds
		mv cell.counts.matrices.rds ${sample_id}.cell.counts.matrices.rds

     	gsutil -q -m cp est_main.log *.rds ${output_directory}/
	}

	output {
		String count_matrices="${output_directory}/${sample_id}.cell.counts.matrices.rds"
		String count_matrix ="${output_directory}/${sample_id}.cell.counts.rds"
		String log="${output_directory}/est_main.log"
	}

	runtime {
		docker: "regevlab/dropest-${dropest_version}"
		preemptible: "${preemptible}"
        zones: zones
		disks: "local-disk " + ceil((3 * size(input_bam, "GB")) + 20)+ " HDD"
		memory :"${memory}"
		cpu:"${cpu}"
	}
}

task run_dropest_umi_correct {
	String memory
	File count_matrix
	Int cpu
	String sample_id
	Int preemptible
	String zones
	String output_directory
	String dropest_version

	command {
		set -e

		monitor_script.sh &

		dropest_umi_correct.R ${count_matrix} ${cpu}

		mv umi_corrected.rds ${sample_id}.umi_corrected.rds
     	gsutil -q -m cp *.rds ${output_directory}/
	}

	output {
		String umi_corrected_count_matrix ="${output_directory}/${sample_id}.umi_corrected.rds"
	}

	runtime {
		docker: "regevlab/dropest-${dropest_version}"
		preemptible: "${preemptible}"
        zones: zones
		disks: "local-disk " + ceil((1.5 * size(count_matrix, "GB")) + 20)+ " HDD"
		memory :"${memory}"
		cpu:"${cpu}"
	}
}








