version 1.0


workflow dropseq_bundle {
  input {
    Int star_cpus = 64
    String star_memory = "57.6G"
    String create_intervals_memory = "3.75G"
    Float star_index_extra_disk_space = 15
    Array[File] fasta_file
    Float fix_gtf_extra_disk_space = 10
    Float add_fasta_prefix_extra_disk_space = 10
    # e.g. --limitGenomeGenerateRAM=124544990592
    String? extra_star_flags
    Array[File] gtf_file
    Int preemptible = 2
    String docker_registry = "cumulusprod"
    Int genomeSAindexNbases = 14
    String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    String drop_seq_tools_version = "2.3.0"
  }


  File? fasta = if (length(fasta_file) > 1) then add_fasta_prefix.fasta else fasta_file[0]
  String docker_registry_stripped = sub(docker_registry, "/+$", "")
  File gtf = fix_gtf.gtf
  String? bundle_name = if (length(fasta_file) > 1) then get_bundle_name.bundle_name else sub(basename(fasta_file[0]), "\.fasta$|\.fa$", "")
  if (length(fasta_file) > 1) {
    call get_bundle_name {
      input:
        input_fasta = fasta_file,
        zones = zones,
        drop_seq_tools_version = drop_seq_tools_version,
        preemptible = preemptible,
        docker_registry = docker_registry_stripped
    }
    call add_fasta_prefix {
      input:
        input_fasta = fasta_file,
        prefix = get_bundle_name.prefix,
        bundle_name = get_bundle_name.bundle_name,
        zones = zones,
        drop_seq_tools_version = drop_seq_tools_version,
        preemptible = preemptible,
        docker_registry = docker_registry_stripped,
        extra_disk_space = add_fasta_prefix_extra_disk_space
    }
  }
  call create_intervals {
    input:
      reduced_gtf = reduce_gtf.reduced_gtf,
      sequence_dictionary = create_sequence_dictionary.dict,
      prefix = bundle_name,
      memory = create_intervals_memory,
      zones = zones,
      drop_seq_tools_version = drop_seq_tools_version,
      preemptible = preemptible,
      docker_registry = docker_registry_stripped
  }
  call create_sequence_dictionary {
    input:
      fasta = fasta,
      output_name = bundle_name + ".dict",
      zones = zones,
      drop_seq_tools_version = drop_seq_tools_version,
      preemptible = preemptible,
      docker_registry = docker_registry_stripped
  }
  call star_index {
    input:
      fasta = fasta,
      gtf = gtf,
      genomeSAindexNbases = genomeSAindexNbases,
      prefix = bundle_name,
      threads = star_cpus,
      memory = star_memory,
      extra_star_flags = extra_star_flags,
      zones = zones,
      drop_seq_tools_version = drop_seq_tools_version,
      preemptible = preemptible,
      docker_registry = docker_registry_stripped,
      extra_disk_space = star_index_extra_disk_space
  }
  call fix_gtf {
    input:
      input_gtf = gtf_file,
      prefix = get_bundle_name.prefix,
      bundle_name = bundle_name,
      preemptible = preemptible,
      zones = zones,
      drop_seq_tools_version = drop_seq_tools_version,
      docker_registry = docker_registry_stripped,
      extra_disk_space = fix_gtf_extra_disk_space
  }
  call reduce_gtf {
    input:
      gtf = gtf,
      dict = create_sequence_dictionary.dict,
      output_name = bundle_name + "_reduced.gtf",
      zones = zones,
      drop_seq_tools_version = drop_seq_tools_version,
      preemptible = preemptible,
      docker_registry = docker_registry_stripped
  }
  call convert_to_ref_flat {
    input:
      gtf = gtf,
      dict = create_sequence_dictionary.dict,
      output_name = bundle_name + ".refFlat",
      zones = zones,
      drop_seq_tools_version = drop_seq_tools_version,
      preemptible = preemptible,
      docker_registry = docker_registry_stripped
  }

  output {
    File ref_flat = convert_to_ref_flat.ref_flat
    File exons_intervals = create_intervals.exons_intervals
    File? output_fasta = if (length(fasta_file) > 1) then add_fasta_prefix.fasta else fasta_file[0]
    File output_gtf = fix_gtf.gtf
    File intergenic_intervals = create_intervals.intergenic_intervals
    File index_tar_gz = star_index.index_tar_gz
    File consensus_intervals = create_intervals.consensus_intervals
    File rrna_intervals = create_intervals.rrna_intervals
    File dict = create_sequence_dictionary.dict
    File genes_intervals = create_intervals.genes_intervals
  }
}
task reduce_gtf {
  input {
    File gtf
    File dict
    String output_name
    String zones
    String drop_seq_tools_version
    Int preemptible
    String docker_registry
  }


  output {
    File reduced_gtf = "${output_name}"
  }
  command <<<

		set -e
		java -Xmx3500m -jar /software/Drop-seq_tools/jar/dropseq.jar ReduceGtf GTF=~{gtf} SD=~{dict} O=~{output_name}

  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + 1 + ceil(size(gtf, "GB") * 2 + size(dict, "GB")) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "3.75 GB"
  }

}
task get_bundle_name {
  input {
    Array[String] input_fasta
    String zones
    String drop_seq_tools_version
    Int preemptible
    String docker_registry
  }


  output {
    String bundle_name = read_string("bundle_name.txt")
    String prefix = read_string("prefix.txt")
  }
  command <<<

		set -e

		python <<CODE
		import os
		from subprocess import check_call
		fasta_files = "~{sep=","  input_fasta}".split(',')
		prefix_list = []
		for i in range(len(fasta_files)):
			prefix = os.path.basename(fasta_files[i])
			dot_index = prefix.rfind('.')
			prefix = prefix[0:dot_index].replace(' ', '_')
			prefix_list.append(prefix)
		with open("bundle_name.txt", "w") as b, open("prefix.txt", "w") as p:
			b.write('_'.join(prefix_list) + '\n')
			p.write(','.join(prefix_list) + '\n')
		CODE


  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk 2 HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "1 GB"
  }

}
task create_sequence_dictionary {
  input {
    File? fasta
    String output_name
    String zones
    String drop_seq_tools_version
    Int preemptible
    String docker_registry
  }


  output {
    File dict = "${output_name}"
  }
  command <<<

		set -e
		java -Xmx3500m -jar /software/picard.jar CreateSequenceDictionary R=~{fasta} O=~{output_name}

  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + 1 + ceil(size(fasta, "GB") * 2) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "3.75G"
  }

}
task convert_to_ref_flat {
  input {
    File gtf
    File dict
    String output_name
    String zones
    String drop_seq_tools_version
    Int preemptible
    String docker_registry
  }


  output {
    File ref_flat = "${output_name}"
  }
  command <<<

		set -e
		java -Xmx3500m -jar /software/Drop-seq_tools/jar/dropseq.jar ConvertToRefFlat ANNOTATIONS_FILE=~{gtf} SD=~{dict} O=~{output_name}

  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + 1 + ceil(size(gtf, "GB") * 2 + size(dict, "GB")) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "3.75 GB"
  }

}
task create_intervals {
  input {
    File reduced_gtf
    String memory
    File sequence_dictionary
    String? prefix
    String zones
    String drop_seq_tools_version
    Int preemptible
    String docker_registry
  }


  output {
    File genes_intervals = "${prefix}.genes.intervals"
    File exons_intervals = "${prefix}.exons.intervals"
    File rrna_intervals = "${prefix}.rRNA.intervals"
    File consensus_intervals = "${prefix}.consensus_introns.intervals"
    File intergenic_intervals = "${prefix}.intergenic.intervals"
  }
  command <<<

		set -e
		mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print $2 }')
		mem="$(($mem/1000))"
		java -Xmx$(echo $mem)m -jar /software/Drop-seq_tools/jar/dropseq.jar CreateIntervalsFiles SD=~{sequence_dictionary} REDUCED_GTF=~{reduced_gtf} OUTPUT=. PREFIX=~{prefix}

  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + 1 + ceil(size(reduced_gtf, "GB") * 2 + size(sequence_dictionary, "GB")) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "${memory}"
  }

}
task add_fasta_prefix {
  input {
    Array[File] input_fasta
    String prefix
    String bundle_name
    String zones
    String drop_seq_tools_version
    Int preemptible
    String docker_registry
    Float extra_disk_space
  }


  output {
    File fasta = "${bundle_name}.fasta"
  }
  command <<<

		set -e

		add_fasta_prefix.py \
		--prefix ~{prefix} \
		--output ~{bundle_name} \
		~{sep=" "  input_fasta}

  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil(extra_disk_space + size(input_fasta[0], "GB") * (length(input_fasta) + 1)) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "3.75 GB"
  }

}
task fix_gtf {
  input {
    Array[File] input_gtf
    String? prefix
    String? bundle_name
    Int preemptible
    String zones
    String drop_seq_tools_version
    String docker_registry
    Float extra_disk_space
  }


  output {
    File gtf = "${bundle_name}.gtf"
  }
  command <<<

		set -e

		fix_gtf.py \
		~{"--prefix " + prefix} \
		--output "~{bundle_name}.gtf" \
		~{sep=" "  input_gtf}

  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil(extra_disk_space + size(input_gtf[0], "GB") * (1 + length(input_gtf))) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: 1
    zones: zones
    memory: "3.75 GB"
  }

}
task star_index {
  input {
    File? fasta
    File gtf
    Int genomeSAindexNbases
    String? extra_star_flags
    String? prefix
    Int threads
    String memory
    String zones
    String drop_seq_tools_version
    Int preemptible
    String docker_registry
    Float extra_disk_space
  }


  output {
    File index_tar_gz = "${prefix}.tgz"
  }
  command <<<

		set -e
		mkdir ~{prefix}
		STAR --runThreadN ~{threads} \
		--runMode genomeGenerate \
		--genomeDir ~{prefix} \
		--genomeFastaFiles ~{fasta} \
		--outFileNamePrefix ~{prefix} \
		--sjdbGTFfile ~{gtf} \
		--genomeSAindexNbases ~{genomeSAindexNbases} \
		~{extra_star_flags}
		tar czf "~{prefix}.tgz" ~{prefix}

  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil(extra_disk_space + size(gtf, "GB") + 22 * size(fasta, "GB")) + " HDD"
    docker: "${docker_registry}/dropseq:${drop_seq_tools_version}"
    cpu: "${threads}"
    zones: zones
    memory: "${memory}"
  }

}

