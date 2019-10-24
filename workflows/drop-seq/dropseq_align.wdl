workflow dropseq_align {
   	String sample_id
   	String star_memory = "57.6G"
    Int star_cpus = 64
    File input_bam

    # e.g. --twopassMode Basic
    String? star_flags
    File star_genome_file
	File refflat
	File gene_intervals
	File genome_fasta
	File genome_dict
	String? zones = "us-east1-d us-west1-a us-west1-b"
	Int preemptible = 2
	String output_directory
	String drop_seq_tools_version
	Int add_bam_tags_disk_space_multiplier = 25
	String? merge_bam_alignment_memory="13G"
	Int? sort_bam_max_records_in_ram = 2000000
    String docker_registry
	call STAR {
		input:
			preemptible=preemptible,
			input_bam=input_bam,
			genome_tar=star_genome_file,
			sample_id=sample_id,
			flags=star_flags,
			memory=star_memory,
			cpu=star_cpus,
			zones=zones,
			output_directory=output_directory,
			drop_seq_tools_version=drop_seq_tools_version,
			docker_registry=docker_registry
	}

	call AddTags {
		input:
			memory=merge_bam_alignment_memory,
			sort_bam_max_records_in_ram=sort_bam_max_records_in_ram,
            cpu=1,
			aligned_bam=STAR.bam,
			unaligned_bam=input_bam,
			sample_id=sample_id,
			preemptible=preemptible,
			genome_fasta=genome_fasta,
			genome_dict=genome_dict,
			refflat=refflat,
			gene_intervals=gene_intervals,
			zones=zones,
			output_directory=output_directory,
			drop_seq_tools_version=drop_seq_tools_version,
			disk_space_multiplier=add_bam_tags_disk_space_multiplier,
			docker_registry=docker_registry
	}


	output {
		String aligned_tagged_bam=AddTags.bam
		String aligned_bam=STAR.bam
		String star_log_final = STAR.star_log_final
		String output_sample_id = AddTags.output_sample_id
	}

}



task STAR {
	String memory
	String? flags
	File input_bam
	Int cpu
	String sample_id
	File genome_tar
	Int preemptible
	String zones
	String output_directory
	String drop_seq_tools_version
	String docker_registry


	command {
		set -o pipefail
		set -e

		mkdir -p genome_dir
		tar xf ${genome_tar} -C genome_dir --strip-components 1

		java -Xmx5000m -jar /software/picard.jar SamToFastq \
		VALIDATION_STRINGENCY=SILENT \
		INPUT=${input_bam} \
		FASTQ=/dev/stdout | \
		STAR \
		--genomeDir genome_dir \
		--outStd Log \
		--readFilesIn /dev/stdin \
		--runThreadN ${cpu} \
		--outSAMtype BAM Unsorted \
		--outFileNamePrefix "${sample_id}_" \
		${" " + flags}

     	gsutil -q -m cp ${sample_id}_* ${output_directory}/
	}

	output {

		String bam="${output_directory}/${sample_id}_Aligned.out.bam"
		String star_log_final="${output_directory}/${sample_id}_Log.final.out"
	}

	runtime {
		docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
		preemptible: "${preemptible}"
        zones: zones
		disks: "local-disk " + ceil(size(genome_tar, "GB")*5 + (3.25 * size(input_bam, "GB")) + 1)+ " HDD"
		memory :"${memory}"
		cpu:"${cpu}"
	}
}


task AddTags {
	String memory
	File aligned_bam
	File unaligned_bam
	String sample_id
	Int preemptible
	File genome_fasta
	File genome_dict
	String zones
	File refflat
	File gene_intervals
	String output_directory
	String drop_seq_tools_version
	Int disk_space_multiplier
	Int cpu
	Int sort_bam_max_records_in_ram
    String docker_registry

	command {
		set -o pipefail
    	set -e

		mkfifo pipe1 pipe2 pipe3

		java -Dsamjdk.compression_level=1 -Xms3000m -jar /software/picard.jar SortSam VALIDATION_STRINGENCY=SILENT \
		INPUT=${aligned_bam} \
		OUTPUT=pipe1 \
		MAX_RECORDS_IN_RAM=${sort_bam_max_records_in_ram} \
		SORT_ORDER=queryname | \
		java -Dsamjdk.compression_level=1 -Xms3000m -jar /software/picard.jar MergeBamAlignment VALIDATION_STRINGENCY=SILENT \
		ALIGNED_BAM=pipe1 \
		UNMAPPED_BAM=${unaligned_bam} \
		OUTPUT=pipe2 \
		REFERENCE_SEQUENCE=${genome_fasta} \
		INCLUDE_SECONDARY_ALIGNMENTS=false \
		PAIRED_RUN=false \
		CLIP_ADAPTERS=false | \
		java -Dsamjdk.compression_level=1 -Xms3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagReadWithInterval VALIDATION_STRINGENCY=SILENT \
		I=pipe2 \
		O=pipe3 \
		INTERVALS=${gene_intervals} \
		TAG=XG | \
		java -Dsamjdk.compression_level=2 -Xms3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagReadWithGeneFunction VALIDATION_STRINGENCY=SILENT \
		INPUT=pipe3 \
		O=${sample_id}_aligned_tagged.bam \
		ANNOTATIONS_FILE=${refflat}

		gsutil -q -m cp ${sample_id}_aligned_tagged.bam ${output_directory}/
	}

	output {

		String bam="${output_directory}/${sample_id}_aligned_tagged.bam"
		String output_sample_id = "${sample_id}"
	}

	runtime {
		docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
		disks: "local-disk " + ceil(20 + size(genome_fasta,"GB") + size(gene_intervals,"GB") + size(genome_dict,"GB") + size(refflat,"GB") +  size(unaligned_bam,"GB")*2 + size(aligned_bam,"GB")*disk_space_multiplier) + " HDD"
		memory :"${memory}"
		preemptible: "${preemptible}"
		zones: zones
		cpu:"${cpu}"
	}
}







