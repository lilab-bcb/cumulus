version 1.0

workflow smartseq2 {
    input {
        # 3-4 columns (entity:sample_id, plate, read1, and optionally read2). gs URL
        File input_tsv_file
        # Output directory, gs URL
        String output_directory
        # Reference to align reads against, GRCh38_ens93filt or GRCm38_ens93filt
        String reference
        # Align reads with 'aligner': hisat2-hca, star, bowtie2 (default: hisat2-hca)
        String aligner = "hisat2-hca"

        # Convert transcript BAM file into genome BAM file
        Boolean output_genome_bam = false
        # smartseq2 version, default to "1.3.0"
        String smartseq2_version = "1.3.0"
        # Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Number of cpus per job
        Int num_cpu = 4
        # Memory string
        String memory = (if aligner != "star" then "3.60G" else "32G")
        # factor to multiply size of R1 and R2 by for RSEM
        Float disk_space_multiplier = 11
        # Number of preemptible tries 
        Int preemptible = 2
        # Disk space for count matrix generation task
        Int generate_count_matrix_disk_space = 10
        # Which docker registry to use: cumulusprod (default) or quay.io/cumulus
        String docker_registry = "cumulusprod"  
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "/+$", "")
    
    Array[Array[String]] data_table = read_tsv(input_tsv_file)


    File acronym_file = "gs://regev-lab/resources/smartseq2/index.tsv"
    # File acronym_file = "smartseq2_index.tsv"
    Map[String, String] acronym2gsurl = read_map(acronym_file)
    # If reference is a url
    Boolean is_url = sub(reference, "^.+\\.(tgz|gz)$", "URL") == "URL"

    String key = reference + "_" + aligner
    File reference_file = (if is_url then reference else acronym2gsurl[key])


    
    scatter (i in range(length(data_table))) {
        Boolean is_se = length(data_table[0]) == 3
        Int pos = i + 1

        # Single-end data
        if (is_se) {
            call run_rsem as run_rsem_se {
                input:
                    reference = reference_file,
                    read1 = data_table[pos][2],
                    sample_name = data_table[pos][0] + "." + data_table[pos][1],
                    aligner = aligner,
                    output_genome_bam = output_genome_bam,
                    smartseq2_version = smartseq2_version,
                    output_directory = output_directory,
                    zones = zones,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space_multiplier = disk_space_multiplier,
                    preemptible = preemptible,
                    docker_registry = docker_registry
            }        
        }
        
        # Paired-end data
        if (!is_se) {
            call run_rsem as run_rsem_pe {
                input:
                    reference = reference_file,
                    read1 = data_table[pos][2],
                    read2 = data_table[pos][3],
                    sample_name = data_table[pos][0] + "." + data_table[pos][1],
                    aligner = aligner,
                    output_genome_bam = output_genome_bam,
                    smartseq2_version = smartseq2_version,
                    output_directory = output_directory,
                    zones = zones,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space_multiplier = disk_space_multiplier,
                    preemptible = preemptible,
                    docker_registry = docker_registry
            }
        }
    }

    call generate_count_matrix {
        input:
            gene_results = select_first([run_rsem_se.rsem_gene, run_rsem_pe.rsem_gene]),
            count_results = select_first([run_rsem_se.rsem_cnt, run_rsem_pe.rsem_cnt]),
            output_directory = output_directory,
            smartseq2_version = smartseq2_version,
            zones = zones,
            memory = memory,
            disk_space = generate_count_matrix_disk_space,
            preemptible = preemptible,
            docker_registry = docker_registry
    }

    output {
        Array[File?] rsem_gene = select_first([run_rsem_se.rsem_gene, run_rsem_pe.rsem_gene])
        Array[File?] rsem_isoform = select_first([run_rsem_se.rsem_isoform, run_rsem_pe.rsem_isoform])
        Array[String?] rsem_trans_bam = select_first([run_rsem_se.rsem_trans_bam, run_rsem_pe.rsem_trans_bam])
        Array[String?] rsem_genome_bam = select_first([run_rsem_se.rsem_genome_bam, run_rsem_pe.rsem_genome_bam])
        Array[File?] rsem_time = select_first([run_rsem_se.rsem_time, run_rsem_pe.rsem_time])
        Array[File?] aligner_log = select_first([run_rsem_se.aligner_log, run_rsem_pe.aligner_log])
        Array[File?] rsem_cnt = select_first([run_rsem_se.rsem_cnt, run_rsem_pe.rsem_cnt])
        Array[File?] rsem_model = select_first([run_rsem_se.rsem_model, run_rsem_pe.rsem_model])
        Array[File?] rsem_theta = select_first([run_rsem_se.rsem_theta, run_rsem_pe.rsem_theta])
        String output_count_matrix = generate_count_matrix.output_count_matrix
    }
}


task run_rsem {
    input {
        File reference
        File read1
        File? read2
        Boolean output_genome_bam
        String sample_name
        String aligner
        String smartseq2_version
        String output_directory
        String zones
        Int num_cpu
        String memory
        Float disk_space_multiplier
        Int preemptible
        String docker_registry  
    }

    Boolean is_star = aligner == "star"
    Boolean is_gzipped = sub(read1, "^.+\\.(gz)$", "GZ") == "GZ"
    Boolean star_gzipped_read_file = is_star && is_gzipped

    command {
        set -e
        export TMPDIR=/tmp

        mkdir -p rsem_ref
        tar xf ~{reference} -C rsem_ref --strip-components 1
        REFERENCE_NAME="$(basename `ls rsem_ref/*.grp` .grp)"
        rsem-calculate-expression --~{aligner} ~{true="--output-genome-bam" false="" output_genome_bam} ~{true="--star-gzipped-read-file" false="" star_gzipped_read_file} ~{true="--paired-end" false="" defined(read2)} -p ~{num_cpu} --append-names --time ~{read1} ~{default="" read2} rsem_ref/$REFERENCE_NAME ~{sample_name}

        gsutil cp ~{sample_name}.transcript.bam ~{output_directory}/
        # mkdir -p ~{output_directory}
        # cp ~{sample_name}.transcript.bam ~{output_directory}/

        if [ -f ~{sample_name}.genome.bam ]
        then
            gsutil cp ~{sample_name}.genome.bam ~{output_directory}/
            # cp ~{sample_name}.genome.bam ~{output_directory}/
        fi
    }

    output {
        File rsem_gene = "~{sample_name}.genes.results"
        File rsem_isoform = "~{sample_name}.isoforms.results"
        String rsem_trans_bam = "~{output_directory}/~{sample_name}.transcript.bam"
        String rsem_genome_bam = if output_genome_bam then "~{output_directory}/~{sample_name}.genome.bam" else ""
        File rsem_time = "~{sample_name}.time"
        File aligner_log = "~{sample_name}.log"
        File rsem_cnt = "~{sample_name}.stat/~{sample_name}.cnt"
        File rsem_model = "~{sample_name}.stat/~{sample_name}.model"
        File rsem_theta = "~{sample_name}.stat/~{sample_name}.theta"
    }

    runtime {
        docker: "~{docker_registry}/smartseq2:~{smartseq2_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(reference, "GB")*5 + (disk_space_multiplier * (size(read1, "GB") + size(read2, "GB"))) + 1)+ " HDD"
        cpu: num_cpu
        preemptible: preemptible
    }
}

task generate_count_matrix {
    input {
        Array[File?] gene_results
        Array[File?] count_results
        String output_directory
        String smartseq2_version
        String zones
        String memory
        Int disk_space
        Int preemptible
        String docker_registry
    }

    command {
        set -e
        export TMPDIR=/tmp

        generate_matrix_ss2.py ~{sep=',' gene_results} count_matrix
        gsutil -m cp -r count_matrix ~{output_directory}/
        # mkdir -p ~{output_directory}
        # cp -r count_matrix ~{output_directory}
    }

    output {
        String output_count_matrix = "~{output_directory}/count_matrix"
    }

    runtime {
        docker: "~{docker_registry}/smartseq2:~{smartseq2_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
    }
}
