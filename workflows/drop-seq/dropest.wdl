workflow dropest {
   	String sample_id
    File input_bam
	String? zones = "us-east1-d us-west1-a us-west1-b"
	Int? preemptible = 2
	String output_directory
	String? dropest_version = "0.8.5"

	String? dropest_memory = "104G"
	Int dropest_cpu = 1
   	File? cellular_barcode_whitelist
	Boolean? apply_directional_umi_correction = true
	# use precise merge strategy (can be slow), recommended to use when the list of real barcodes is not available
	Boolean? merge_barcodes_precise = true

	# Minimal number of genes for cells after the merge procedure
	Int? genes_min = 100
	# maximal number of output cells
	Int? cells_max

	# save separate count matrices for exons, introns and exon/intron spanning reads
	Boolean? velocyto = true

	#  Threshold for the merge procedure
	Float min_merge_fraction  = 0.2
	# Max edit distance between barcodes.
	Int? max_cb_merge_edit_distance = 2

	# Max edit distance between UMIs
	Int? max_umi_merge_edit_distance = 1

	# Minimal number of genes for cells before the merge procedure. Used mostly for optimization.
	Int? min_genes_before_merge = 20


	String? cellular_barcode_tag = "XC"
	String? umi_tag = "XM"
    String? gene_id_tag = "XG"
    String? cellular_barcode_quality_tag = "CY"
    String? umi_quality_tag = "UY"
    String docker_registry

	call run_dropest {
		input:
			memory=dropest_memory,
			cpu=dropest_cpu,
			sample_id=sample_id,
			input_bam=input_bam,
			genes_min=genes_min,
			cells_max=cells_max,
			min_merge_fraction=min_merge_fraction,
			max_cb_merge_edit_distance=max_cb_merge_edit_distance,
			max_umi_merge_edit_distance=max_umi_merge_edit_distance,
			min_genes_before_merge=min_genes_before_merge,
			velocyto=velocyto,
			cellular_barcode_whitelist=cellular_barcode_whitelist,
			apply_directional_umi_correction=apply_directional_umi_correction,
			merge_barcodes_precise=merge_barcodes_precise,
			cellular_barcode_tag =cellular_barcode_tag,
			umi_tag = umi_tag,
			gene_id_tag = gene_id_tag,
			cellular_barcode_quality_tag = cellular_barcode_quality_tag,
			umi_quality_tag = umi_quality_tag,
			preemptible=preemptible,
			zones=zones,
			output_directory=output_directory,
			dropest_version=dropest_version,
			docker_registry=docker_registry
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
		String? count_matrices = run_dropest.count_matrices
		String count_matrix = run_dropest.count_matrix
		String bam = run_dropest.bam
		String log= run_dropest.log
#		String? umi_corrected_count_matrix = run_dropest_umi_correct.umi_corrected_count_matrix

	}

}



task run_dropest {
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
	command {
		set -e

		create_dropest_config.py \
		--output dropest_config.xml \
		--min_merge_fraction ${min_merge_fraction} \
        --max_cb_merge_edit_distance ${max_cb_merge_edit_distance} \
        --max_umi_merge_edit_distance ${max_umi_merge_edit_distance} \
        --min_genes_before_merge ${min_genes_before_merge} \
        --cb ${cellular_barcode_tag} \
        --umi ${umi_tag} \
        --gene ${gene_id_tag} \
        --cb_quality ${cellular_barcode_quality_tag} \
        --umi_quality ${umi_quality_tag} \
        ${"--whitelist " + cellular_barcode_whitelist}

		dropest -F \
		${true="-V " false='' velocyto} \
		${true="-u " false='' apply_directional_umi_correction} \
		${true="-M " false='' merge_barcodes_precise} \
		${"-C " + cells_max} \
		${"-G " + genes_min} \
		-c dropest_config.xml \
		-f ${input_bam}

		mv ${input_bam_name}.filtered.bam ${input_bam_name}.dropest_filtered.bam
		mv cell.counts.rds ${sample_id}.cell.counts.rds


		if [ '${velocyto}' == 'true' ]; then
			mv cell.counts.matrices.rds ${sample_id}.cell.counts.matrices.rds
		fi

     	gsutil -q -m cp est_main.log *.rds *.bam ${output_directory}/
	}

	output {
		String? count_matrices="${output_directory}/${sample_id}.cell.counts.matrices.rds"
		String count_matrix ="${output_directory}/${sample_id}.cell.counts.rds"
		String bam = "${output_directory}/${input_bam_name}.dropest_filtered.bam"
		String log="${output_directory}/est_main.log"
	}

	runtime {
		docker: "${docker_registry}/dropest:${dropest_version}"
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
    String docker_registry

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
		docker: "${docker_registry}/dropest:${dropest_version}"
		preemptible: "${preemptible}"
        zones: zones
		disks: "local-disk " + ceil((1.5 * size(count_matrix, "GB")) + 20)+ " HDD"
		memory :"${memory}"
		cpu:"${cpu}"
	}
}








