workflow dropseq_bundle {

#   Reference genome fasta file(s)
    Array[File] fasta_file

#   Reference genome gtf file(s)
#   If more than one species, fasta and gtf files need to be in the same order
    Array[File] gtf_file

#   length (bases) of the SA pre-indexing string. Typically between 10 and 15.
#   Longer strings will use much more memory, but allow faster searches.
#   For small genomes, the parameter --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 1)
    Int genomeSAindexNbases = 14
	Int star_cpus = 64
	String star_memory = "57.6G"
    String? multi_species_bundle_name
	String bundle_name = if(length(fasta_file)>1) then multi_species_bundle_name else sub(basename(fasta_file[0]), "\\.fasta$|\\.fa$", "")
	File gtf = fix_gtf.fixed_gtf
	File? fasta = if(length(fasta_file)>1) then add_fasta_prefix.merged_fasta else fasta_file[0]
	String? zones = "us-east1-d us-west1-a us-west1-b"
	# Number of preemptible tries
    Int? preemptible = 2
    if(length(fasta_file)>1) {
        call add_fasta_prefix {
            input:
                fasta=fasta_file,
                output_name=bundle_name + '.fasta'
        }
    }

    call fix_gtf {
        input:
            gtf=gtf_file,
            output_name=bundle_name + '.gtf',
            fasta=fasta_file
    }


    call create_sequence_dictionary {
        input:
            fasta=fasta,
            output_name = bundle_name + ".dict"
    }
    call convert_to_ref_flat {
        input:
            gtf=gtf,
            dict=create_sequence_dictionary.dict,
            output_name=bundle_name +  ".refFlat"
    }
    call reduce_gtf {
        input:
            gtf=gtf,
            dict=create_sequence_dictionary.dict,
            output_name=bundle_name +  "_reduced.gtf"
    }
    call create_intervals {
        input:
            reduced_gtf=reduce_gtf.reduced_gtf,
            sequence_dictionary=create_sequence_dictionary.dict,
            prefix=bundle_name
    }
    call star_index {
        input:
			fasta=fasta,
			gtf=gtf,
			genomeSAindexNbases=genomeSAindexNbases,
			threads=star_cpus,
			memory=star_memory,
			prefix=bundle_name
    }

    output {
		File genes_intervals =create_intervals.genes_intervals
		File exons_intervals = create_intervals.exons_intervals
		File rrna_intervals = create_intervals.rrna_intervals
		File consensus_intervals = create_intervals.consensus_intervals
		File intergenic_intervals = create_intervals.intergenic_intervals
		File dict = create_sequence_dictionary.dict
		File output_fastq = if(length(fasta_file)>1) then add_fasta_prefix.merged_fasta else fasta_file[0]
		File index_tar_gz = star_index.index_tar_gz
		File ref_flat = convert_to_ref_flat.ref_flat
    }
}

task add_fasta_prefix {
    Array[File] fasta
    String output_name
	Int preemptible
	String zones
	
    command {
        add_fasta_prefix.py --output ${output_name} ${sep=' ' fasta}
    }
    output {
        File merged_fasta="${output_name}"
    }
    runtime {
		docker: "regevlab/dropseq-0.0.1"
		preemptible: "${preemptible}"
		zones: zones
		memory: "2 GB"
		disks: "local-disk " + ceil(size(fasta[0],"GB")*(length(fasta)+1)) + " HDD"
		cpu: 1
    }
}

task fix_gtf {
    Array[File] gtf
    String output_name
    Array[String] fasta_file
    Int preemptible
    String zones
    	
	
    command {
        fix_gtf.py --prefix ${sep=',' fasta_file} --output ${output_name} ${sep=' ' gtf}
    }
    output {
        File fixed_gtf="${output_name}"
    }
    runtime {
		docker: "regevlab/dropseq-0.0.1"
		preemptible: "${preemptible}"
		zones: zones
		memory: "2 GB"
		disks: "local-disk " + ceil(size(gtf[0],"GB")*(1+length(gtf))) + " HDD"
		cpu: 1
  }
}

task create_sequence_dictionary {
    File fasta
    String output_name
	Int preemptible
	String zones
	
    command {
        java -Xmx2g -jar /home/picard-tools-1.141/picard.jar CreateSequenceDictionary R=${fasta} O=${output_name}
    }

    output {
        File dict="${output_name}"
    }
    runtime {
		docker: "regevlab/dropseq-0.0.1"
		preemptible: "${preemptible}"
		zones: zones
		memory: "2 GB"
		disks: "local-disk " + ceil(size(fasta,"GB")*2)  + " HDD"
		cpu: 1
  }
}

task convert_to_ref_flat {
    File gtf
    File dict
    String output_name
	Int preemptible
	String zones
	
    command {
        java -Xmx4g -jar /home/Drop-seq_tools-1.13/jar/dropseq.jar ConvertToRefFlat ANNOTATIONS_FILE=${gtf} SEQUENCE_DICTIONARY=${dict} O=${output_name}
    }
    output {
        File ref_flat="${output_name}"
    }
    runtime {
		docker: "regevlab/dropseq-0.0.1"
		preemptible: "${preemptible}"
		zones: zones
		memory: "4 GB"
		disks: "local-disk " + ceil(size(gtf,"GB")*2 + size(dict,"GB"))  + " HDD"
		cpu: 1
    }
}

task reduce_gtf {
    File gtf
    File dict
    String output_name
	Int preemptible
	String zones
	
    command {
        java -Xmx4g -jar /home/Drop-seq_tools-1.13/jar/dropseq.jar ReduceGTF GTF=${gtf} SEQUENCE_DICTIONARY=${dict} O=${output_name}
    }
    output {
        File reduced_gtf="${output_name}"
    }
    runtime {
		docker: "regevlab/dropseq-0.0.1"
		preemptible: "${preemptible}"
		zones: zones
		memory: "4 GB"
		disks: "local-disk " + ceil(size(gtf,"GB")*2 + size(dict,"GB"))  + " HDD"
		cpu: 1
    }
}

task create_intervals {
    File reduced_gtf
    File sequence_dictionary
    String prefix
	Int preemptible
	String zones
	
    command {
        java -Xmx4g -jar /home/Drop-seq_tools-1.13/jar/dropseq.jar CreateIntervalsFiles SD=${sequence_dictionary} REDUCED_GTF=${reduced_gtf} OUTPUT=. PREFIX=${prefix}
    }
    output {
        File genes_intervals = "${prefix}.genes.intervals"
        File exons_intervals = "${prefix}.exons.intervals"
        File rrna_intervals = "${prefix}.rRNA.intervals"
        File consensus_intervals = "${prefix}.consensus_introns.intervals"
        File intergenic_intervals = "${prefix}.intergenic.intervals"
    }

    runtime {
		docker: "regevlab/dropseq-0.0.1"
		preemptible: "${preemptible}"
		zones: zones
		memory: "4 GB"
		disks: "local-disk " + ceil(size(reduced_gtf,"GB")*2 + size(sequence_dictionary,"GB"))  + " HDD"
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
	Int preemptible
	String zones

    command {
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
        File index_tar_gz ="${prefix}.tar.gz"
    }
    runtime {
       	docker: "regevlab/dropseq-0.0.1"
		preemptible: "${preemptible}"
		zones: zones
        memory: "${memory}"
        disks: "local-disk " + ceil(5 + size(gtf,"GB") + 20*size(fasta,"GB"))  + " HDD"
        cpu: "${threads}"
  }
}



