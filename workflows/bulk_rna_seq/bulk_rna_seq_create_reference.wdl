version 1.0

workflow bulk_rna_seq_create_reference {
    input {
        # Input Genome FASTA file
        File fasta
        # Input Gene annotation file in GTF format
        File gtf
        # Output reference name 
        String reference_name
        Int cpu = 32
        String memory = "40G"
        Int extra_disk_space = 120
        # Number of preemptible tries 
        Int preemptible = 2
        String rsem_docker = "cumulusprod/bulk_rna_seq:1.0.0"
    }


    call rsem_prepare_reference {
        input:
            fasta=fasta,
            gtf=gtf,
            reference_name = reference_name,
            aligner = "star",
            preemptible=preemptible,
            cpu=cpu,
            memory=memory,
            extra_disk_space=extra_disk_space,
            docker=rsem_docker
    }

    output {
        File output_reference = rsem_prepare_reference.output_reference
    }
}

task rsem_prepare_reference {
    input {
        File fasta
        File gtf

        String reference_name
        String aligner
        Int preemptible
        Int cpu
        String memory
        Int extra_disk_space
        String docker
    }

    command {
        set -e

        monitor_script.sh > monitoring.log &

        mkdir ~{reference_name}_~{aligner}
        rsem-prepare-reference --gtf ~{gtf} --~{aligner} -p ~{cpu} ~{fasta} ~{reference_name}_~{aligner}/rsem_ref
        tar -czf ~{reference_name}_~{aligner}.tar.gz ~{reference_name}_~{aligner}


    }

    output {
        File output_reference = "~{reference_name}_~{aligner}.tar.gz"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        disks: "local-disk " + ceil(extra_disk_space + size(gtf, "GB") + size(fasta, "GB")) + " HDD"
        docker: docker
        preemptible: preemptible
        cpu: cpu
        memory: memory
    }
}
