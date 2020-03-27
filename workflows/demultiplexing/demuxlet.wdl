version 1.0

workflow demuxlet {
    input {
        String sample_id
        String output_directory
        String input_rna
        String input_bam
        String ref_genotypes
        Int min_num_genes
        Int min_MQ = 20
        String alpha = "0.1,0.2,0.3,0.4,0.5"
        Int min_TD = 0
        String tag_group = "CB"
        String tag_UMI  = "UB"
        String field = "GT"
        Float geno_error = 0.1

        String zones = "us-central1-b us-east1-d us-west1-a us-west1-b"
        Int num_cpu = 1
        Int memory = 10
        Int extra_disk_space = 2
        Int preemptible = 2
        String docker_registry = "cumulusprod"
        # gs://fc-xx/out/filtered_feature_bc_matrix.h5
        # gs://fc-xx/out/possorted_genome_bam.bam
        # gs://fc-xx/out/possorted_genome_bam.bam.bai
        # gs://fc-xx/out/filtered_feature_bc_matrix/barcodes.tsv.gz
        
    }

    call demuxlet_task {
        input:
            sample_id = sample_id,
            output_directory = output_directory,
            input_rna = input_rna,
            input_bam = input_bam,
            ref_genotypes = ref_genotypes,
            min_num_genes = min_num_genes,
            min_MQ = min_MQ,
            alpha = alpha,
            min_TD = min_TD,
            tag_group = tag_group,
            tag_UMI = tag_UMI,
            field = field,
            geno_error = geno_error,
            docker_registry = docker_registry,
            num_cpu = num_cpu,
            memory = memory,
            extra_disk_space = extra_disk_space,
            zones = zones,
            preemptible = preemptible
    }

    output {
        File output_tsv  = demuxlet_task.output_tsv
        File monitoringLog = demuxlet_task.monitoringLog
    }

}

task demuxlet_task {
    input {
        String sample_id
        String output_directory
        File input_rna
        File input_bam
        File ref_genotypes
        Int min_num_genes
        Int min_MQ
        String alpha
        Int min_TD
        String tag_group
        String tag_UMI
        String field
        Float geno_error

        String docker_registry
        Int num_cpu
        Int memory
        Int extra_disk_space
        Int preemptible
        String zones
    }

    command {
        set -e
        monitor_script.sh > monitoring.log &

        python /software/extract_barcodes_from_rna.py ~{input_rna} ~{sample_id}.barcodes.tsv ~{min_num_genes}

        python <<CODE
        from subprocess import check_call

        call_args = ['popscle', 'demuxlet', '--sam', '~{input_bam}', '--vcf', '~{ref_genotypes}', '--group-list', '~{sample_id}.barcodes.tsv', '--field', '~{field}', '--out', '~{sample_id}']

        if '~{min_MQ}' is not '':
            call_args.extend(['--min-MQ', '~{min_MQ}'])
        if '~{min_TD}' is not '':
            call_args.extend(['--min-TD', '~{min_TD}'])
        if '~{alpha}' is not '':
            alpha_list = '~{alpha}'.split(',')
            prefix_list = ['--alpha'] * len(alpha_list)
            alpha_args = list(sum(list(zip(prefix_list, alpha_list)), ()))
            call_args.extend(alpha_args)
        if '~{tag_group}' is not '':
            call_args.extend(['--tag-group', '~{tag_group}'])
        if '~{tag_UMI}' is not '':
            call_args.extend(['--tag-UMI', '~{tag_UMI}'])
        if '~{geno_error}' is not '':
            call_args.extend(['--geno-error-offset', '~{geno_error}'])

        print(' '.join(call_args))
        check_call(call_args)
        CODE

    }

    output {
        File output_tsv  = "~{sample_id}.best"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/demuxlet:0.1-beta"
        zones: zones
        memory: "~{memory}G"
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(extra_disk_space + size(input_bam,"GB") +  size(ref_genotypes,"GB"))  + " HDD"
        cpu: num_cpu
        preemptible: preemptible
    }
}

