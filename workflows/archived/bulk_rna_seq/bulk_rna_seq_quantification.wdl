version 1.0

workflow bulk_rna_seq_quantification {
    input {
        Array[File] read1
        Array[File] read2
        String sample_name
        String reference
        Boolean output_genome_bam = false
        Int num_cpu = 4
        String memory = "32G"
        String aligner = "star"
        # factor to multiply size of R1 and R2 by
        Float disk_space_multiplier = 4
        Float extra_disk_space = 2
        Int preemptible = 2
        String docker = "cumulusprod/smartseq2:1.1.0"
    }

    File acronym_file = "gs://regev-lab/resources/smartseq2/index.tsv"
    Map[String, String] acronym2gsurl = read_map(acronym_file)
    # If reference is a url
    Boolean is_url = sub(reference, "^.+\\.(tgz|gz)$", "URL") == "URL"

    String key = reference + "_" + aligner
    File reference_file = (if is_url then reference else acronym2gsurl[key])
    call run_rsem {
        input:
            reference = reference_file,
            read1 = read1,
            read2 = read2,
            extra_disk_space = extra_disk_space,
            sample_name = sample_name,
            aligner = aligner,
            output_genome_bam = output_genome_bam,
            num_cpu = num_cpu,
            memory = memory,
            disk_space_multiplier = disk_space_multiplier,
            preemptible = preemptible,
            docker = docker
    }

    output {
        File rsem_gene = run_rsem.rsem_gene
        File rsem_isoform = run_rsem.rsem_isoform
        File rsem_trans_bam = run_rsem.rsem_trans_bam
        File rsem_time = run_rsem.rsem_time
        File aligner_log = run_rsem.aligner_log
        File rsem_cnt = run_rsem.rsem_cnt
        File rsem_model = run_rsem.rsem_model
        File rsem_theta = run_rsem.rsem_theta
        Array[File] rsem_genome_bam = run_rsem.rsem_genome_bam
    }
}

task run_rsem {
    input {
        File reference
        Array[File] read1
        Array[File] read2
        Boolean output_genome_bam
        String sample_name
        String aligner
        Float extra_disk_space
        Int num_cpu
        String memory
        Float disk_space_multiplier
        Int preemptible
        String docker
    }

    Boolean is_star = aligner == "star"
    Boolean paired_end = length(read2) > 0
    Boolean is_gzipped = sub(read1[0], "^.+\\.(gz)$", "GZ") == "GZ"
    Boolean star_gzipped_read_file = is_star && is_gzipped

    command {
        set -e

        monitor_script.sh &

        mkdir -p rsem_ref
        tar xf ~{reference} -C rsem_ref --strip-components 1
        REFERENCE_NAME="$(basename `ls rsem_ref/*.grp` .grp)"
        rsem-calculate-expression --~{aligner} ~{true="--output-genome-bam" false="" output_genome_bam} ~{true="--star-gzipped-read-file" false="" star_gzipped_read_file} ~{true="--paired-end" false="" paired_end} -p ~{num_cpu} --append-names --time ~{sep=',' read1} ~{sep=',' read2} rsem_ref/$REFERENCE_NAME ~{sample_name}
    }

    output {
        File rsem_gene = "~{sample_name}.genes.results"
        File rsem_isoform = "~{sample_name}.isoforms.results"
        File rsem_trans_bam = "~{sample_name}.transcript.bam"
        File rsem_time = "~{sample_name}.time"
        File aligner_log = "~{sample_name}.log"
        File rsem_cnt = "~{sample_name}.stat/~{sample_name}.cnt"
        File rsem_model = "~{sample_name}.stat/~{sample_name}.model"
        File rsem_theta = "~{sample_name}.stat/~{sample_name}.theta"
        Array[File] rsem_genome_bam = glob("~{sample_name}.genome.bam")
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(reference, "GB")*5 + (disk_space_multiplier * (size(read1, "GB") + size(read2, "GB"))) + extra_disk_space)+ " HDD"
        cpu: num_cpu
        preemptible: preemptible
    }
}


