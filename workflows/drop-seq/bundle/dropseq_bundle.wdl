workflow dropseq_bundle {
#	Reference genome fasta file(s)
	Array[File] fasta_file

#	Reference genome gtf file(s)
#	If more than one species, fasta and gtf files need to be in the same order
	Array[File] gtf_file

#	length (bases) of the SA pre-indexing string. Typically between 10 and 15.
#	Longer strings will use much more memory, but allow faster searches.
#	For small genomes, the parameter --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 1)
	Int genomeSAindexNbases = 14
	Int? star_cpus = 64
	String? star_memory = "57.6G"
	String? zones = "us-east1-d us-west1-a us-west1-b"
	String? drop_seq_tools_version = "2.3.0"
	Int? preemptible = 2
	if(length(fasta_file)>1) {
		call get_bundle_name {
			input:
				input_fasta=fasta_file,
				drop_seq_tools_version=drop_seq_tools_version,
				zones=zones,
				preemptible=preemptible
		}

		call add_fasta_prefix {
			input:
				input_fasta=fasta_file,
				prefix=get_bundle_name.prefix,
				bundle_name=get_bundle_name.bundle_name,
				drop_seq_tools_version=drop_seq_tools_version,
				zones=zones,
				preemptible=preemptible
		}
	}
	String? bundle_name = if(length(fasta_file)>1) then get_bundle_name.bundle_name else sub(basename(fasta_file[0]), "\\.fasta$|\\.fa$", "")
	File? fasta = if(length(fasta_file)>1) then add_fasta_prefix.fasta else fasta_file[0]

	call fix_gtf {
		input:
			input_gtf=gtf_file,
			prefix=get_bundle_name.prefix,
			bundle_name=bundle_name,
			drop_seq_tools_version=drop_seq_tools_version,
			zones=zones,
			preemptible=preemptible
	}

	File gtf = fix_gtf.gtf

	call create_sequence_dictionary {
		input:
			fasta=fasta,
			output_name = bundle_name + ".dict",
			drop_seq_tools_version=drop_seq_tools_version,
			zones=zones,
			preemptible=preemptible
	}
	call convert_to_ref_flat {
		input:
			gtf=gtf,
			dict=create_sequence_dictionary.dict,
			output_name=bundle_name +  ".refFlat",
			drop_seq_tools_version=drop_seq_tools_version,
			zones=zones,
			preemptible=preemptible
	}
	call reduce_gtf {
		input:
			gtf=gtf,
			dict=create_sequence_dictionary.dict,
			output_name=bundle_name +  "_reduced.gtf",
			drop_seq_tools_version=drop_seq_tools_version,
			zones=zones,
			preemptible=preemptible
	}
	call create_intervals {
		input:
			reduced_gtf=reduce_gtf.reduced_gtf,
			sequence_dictionary=create_sequence_dictionary.dict,
			prefix=bundle_name,
			drop_seq_tools_version=drop_seq_tools_version,
			zones=zones,
			preemptible=preemptible
	}
	call star_index {
		input:
			fasta=fasta,
			gtf=gtf,
			genomeSAindexNbases=genomeSAindexNbases,
			threads=star_cpus,
			memory=star_memory,
			prefix=bundle_name,
			drop_seq_tools_version=drop_seq_tools_version,
			zones=zones,
			preemptible=preemptible
	}

	output {
		File genes_intervals = create_intervals.genes_intervals
		File exons_intervals = create_intervals.exons_intervals
		File rrna_intervals = create_intervals.rrna_intervals
		File consensus_intervals = create_intervals.consensus_intervals
		File intergenic_intervals = create_intervals.intergenic_intervals
		File dict = create_sequence_dictionary.dict
		File? output_fasta = if(length(fasta_file)>1) then add_fasta_prefix.fasta else fasta_file[0]
		File output_gtf = fix_gtf.gtf
		File index_tar_gz = star_index.index_tar_gz
		File ref_flat = convert_to_ref_flat.ref_flat
	}
}

task get_bundle_name {
	Array[String] input_fasta
	String zones
	String drop_seq_tools_version
	Int preemptible

	command {
		set -e

		python <<CODE
		import os
		from subprocess import check_call
		fasta_files = "${sep=',' input_fasta}".split(',')
		prefix_list = []
		for i in range(len(fasta_files)):
			prefix = os.path.basename(fasta_files[i])
			dot_index = prefix.rfind('.')
			prefix = prefix[0:dot_index]
			prefix_list.append(prefix)
		with open("bundle_name.txt", "w") as b, open("prefix.txt", "w") as p:
			b.write('\n'.join(prefix_list) + '\n')
			p.write(','.join(prefix_list) + '\n')
		CODE

	}
	output {
		String bundle_name = read_string("bundle_name.txt")
		String prefix = read_string("prefix.txt")

	}
	runtime {
		docker: "cumulusprod/dropseq:${drop_seq_tools_version}"
		preemptible: "${preemptible}"
		zones: zones
		memory: "1 GB"
		disks: "local-disk 2 HDD"
		cpu: 1
	}
}

task add_fasta_prefix {
	Array[File] input_fasta
	String prefix
	String bundle_name
	String zones
	String drop_seq_tools_version
	Int preemptible
	
	command {
		set -e
		add_fasta_prefix.py', '--prefix', ','.join(prefix_list), '--output', '_'.join(prefix_list), ${sep=' ' input_fasta}])


	}
	output {
		File fasta="${bundle_name}.fasta"

	}
	runtime {
		docker: "cumulusprod/dropseq:${drop_seq_tools_version}"
		preemptible: "${preemptible}"
		zones: zones
		memory: "3.75 GB"
		disks: "local-disk " + 10+ceil(size(input_fasta[0],"GB")*(length(input_fasta)+1)) + " HDD"
		cpu: 1
	}
}


task fix_gtf {
	Array[File] input_gtf
	String? prefix
	String bundle_name
	Int preemptible
	String zones
	String drop_seq_tools_version

	command {
		set -e
		fix_gtf.py \
		${"--prefix" + prefix} \
		--output ${bundle_name}.gtf \
		${sep=' ' input_gtf}
	}
	output {
		File gtf="${bundle_name}.gtf"
	}
	runtime {
		docker: "cumulusprod/dropseq:${drop_seq_tools_version}"
		preemptible: "${preemptible}"
		zones: zones
		memory: "3.75 GB"
		disks: "local-disk " + ceil(size(input_gtf[0],"GB")*(1+length(input_gtf))) + " HDD"
		cpu: 1
  }
}

task create_sequence_dictionary {
	File fasta
	String output_name
	String zones
	String drop_seq_tools_version
	Int preemptible
	
	command {
		set -e
		java -Xmx3500m -jar /software/picard.jar CreateSequenceDictionary R=${fasta} O=${output_name}
	}

	output {
		File dict="${output_name}"
	}
	runtime {
		docker: "cumulusprod/dropseq:${drop_seq_tools_version}"
		preemptible: "${preemptible}"
		zones: zones
		memory: "3.75G"
		disks: "local-disk " + 1+ceil(size(fasta,"GB")*2)  + " HDD"
		cpu: 1
	}
}

task convert_to_ref_flat {
	File gtf
	File dict
	String output_name
	String zones
	String drop_seq_tools_version
	Int preemptible
	
	command {
		set -e
		java -Xmx3500m -jar /software/Drop-seq_tools/jar/dropseq.jar ConvertToRefFlat ANNOTATIONS_FILE=${gtf} SD=${dict} O=${output_name}
	}
	output {
		File ref_flat="${output_name}"
	}
	runtime {
		docker: "cumulusprod/dropseq:${drop_seq_tools_version}"
		preemptible: "${preemptible}"
		zones: zones
		memory: "3.75 GB"
		disks: "local-disk " + 1+ceil(size(gtf,"GB")*2 + size(dict,"GB"))  + " HDD"
		cpu: 1
	}
}

task reduce_gtf {
	File gtf
	File dict
	String output_name
	String zones
	String drop_seq_tools_version
	Int preemptible
	
	command {
		set -e
		java -Xmx3500m -jar /software/Drop-seq_tools/jar/dropseq.jar ReduceGtf GTF=${gtf} SD=${dict} O=${output_name}
	}
	output {
		File reduced_gtf="${output_name}"
	}
	runtime {
		docker: "cumulusprod/dropseq:${drop_seq_tools_version}"
		preemptible: "${preemptible}"
		zones: zones
		memory: "3.75 GB"
		disks: "local-disk " + 1+ceil(size(gtf,"GB")*2 + size(dict,"GB"))  + " HDD"
		cpu: 1
	}
}

task create_intervals {
	File reduced_gtf
	File sequence_dictionary
	String prefix
	String zones
	String drop_seq_tools_version
	Int preemptible
	
	command {
		set -e
		java -Xmx3500m -jar /software/Drop-seq_tools/jar/dropseq.jar CreateIntervalsFiles SD=${sequence_dictionary} REDUCED_GTF=${reduced_gtf} OUTPUT=. PREFIX=${prefix}
	}
	output {
		File genes_intervals = "${prefix}.genes.intervals"
		File exons_intervals = "${prefix}.exons.intervals"
		File rrna_intervals = "${prefix}.rRNA.intervals"
		File consensus_intervals = "${prefix}.consensus_introns.intervals"
		File intergenic_intervals = "${prefix}.intergenic.intervals"
	}

	runtime {
		docker: "cumulusprod/dropseq:${drop_seq_tools_version}"
		preemptible: "${preemptible}"
		zones: zones
		memory: "3.75 GB"
		disks: "local-disk " + 1+ceil(size(reduced_gtf,"GB")*2 + size(sequence_dictionary,"GB"))  + " HDD"
		cpu: 1
	}
}

task star_index {
	File fasta
	File gtf
	Int genomeSAindexNbases
	String prefix
	Int threads
	String memory
	String zones
	String drop_seq_tools_version
	Int preemptible

	command {
		set -e
		mkdir ${prefix}
		STAR --runThreadN ${threads} \
		--runMode genomeGenerate \
		--genomeDir ${prefix} \
		--genomeFastaFiles ${fasta} \
		--outFileNamePrefix ${prefix} \
		--sjdbGTFfile ${gtf} \
		--genomeSAindexNbases ${genomeSAindexNbases}
		tar czf "${prefix}.tgz" ${prefix}
	}
	output {
		File index_tar_gz ="${prefix}.tgz"
	}
	runtime {
		docker: "cumulusprod/dropseq:${drop_seq_tools_version}"
		preemptible: "${preemptible}"
		zones: zones
		memory: "${memory}"
		disks: "local-disk " + ceil(5 + size(gtf,"GB") + 20*size(fasta,"GB"))  + " HDD"
		cpu: "${threads}"
  }
}

