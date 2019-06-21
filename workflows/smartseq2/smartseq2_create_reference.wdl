workflow smartseq2_create_reference {
	File fasta
	File gtf
	# smartseq2 version, default to "0.1.0"
	String? smartseq2_version = "0.1.0"
	String? zones = "us-central1-b"
	Int? cpu = 8
	String? memory = "7.2G"
	# disk space in GB
	Int? extra_disk_space = 15
	# Number of preemptible tries 
	Int? preemptible = 2

	call rsem_prepare_reference {
		input:
			fasta=fasta,
			gtf=gtf,
			smartseq2_version=smartseq2_version,
			zones=zones,
			preemptible=preemptible,
			cpu=cpu,
			memory=memory,
			disk_space=extra_disk_space
	}
}

task rsem_prepare_reference {
	File fasta
    File gtf
	String smartseq2_version
	String zones
	Int preemptible
	Int cpu
	String memory
	Int disk_space

	command {
		set -e

		mkdir rsem_ref
		rsem-prepare-reference --gtf ${gtf} --bowtie2 -p ${cpu} ${fasta} rsem_ref/rsem_ref
		tar -czf rsem_ref.tar.gz rsem_ref
	}

	output {
		File reference = "rsem_ref.tar.gz"
	}

	runtime {
		disks: "local-disk " + ceil(disk_space + 8*size(fasta,"GB") + size(gtf,"GB")) + " HDD"
		docker: "regevlab/smartseq2-${smartseq2_version}"
		zones: zones
		preemptible: "${preemptible}"
		cpu:"${cpu}"
		memory:"${memory}"
	}
}

