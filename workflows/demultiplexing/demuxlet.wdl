version 1.0

workflow demuxlet {
    input {
        String sample_id
        String output_directory
        String input_rna
        String input_bam
        String ref_genotypes
        Int min_num_genes
        String field = "GT"
        Int? min_MQ
        String? alpha
        Int? min_TD
        String? tag_group
        String? tag_UMI
        Float? geno_error

        String zones = "us-central1-b us-east1-d us-west1-a us-west1-b"
        Int num_cpu = 1
        Int memory = 10
        Int extra_disk_space = 2
        Int preemptible = 2
        String docker_registry = "cumulusprod"
        String demuxlet_version = "0.1b"
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
            demuxlet_version = demuxlet_version,
            num_cpu = num_cpu,
            memory = memory,
            extra_disk_space = extra_disk_space,
            zones = zones,
            preemptible = preemptible
    }

    output {
        String output_folder = demuxlet_task.output_folder
        File output_zarr = demuxlet_task.output_zarr
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
        String field
        Int? min_MQ
        String? alpha
        Int? min_TD
        String? tag_group
        String? tag_UMI
        Float? geno_error

        String docker_registry
        String demuxlet_version
        Int num_cpu
        Int memory
        Int extra_disk_space
        Int preemptible
        String zones
    }

    command {
        set -e
        monitor_script.sh > monitoring.log &

        mkdir result
        python /software/extract_barcodes_from_rna.py ~{input_rna} ~{sample_id}.barcodes.tsv ~{min_num_genes}

        python <<CODE
        from subprocess import check_call

        call_args = ['popscle', 'demuxlet', '--sam', '~{input_bam}', '--vcf', '~{ref_genotypes}', '--group-list', '~{sample_id}.barcodes.tsv', '--field', '~{field}', '--out', 'result/~{sample_id}']

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

        python /software/generate_zarr.py result/~{sample_id}.best ~{input_rna} result/~{sample_id}_demux.zarr.zip

        gsutil -q -m rsync -r result ~{output_directory}/~{sample_id}
        # mkdir -p ~{output_directory}/~{sample_id}
        # cp result/* ~{output_directory}/~{sample_id}
    }

    output {
        String output_folder = "~{output_directory}/~{sample_id}"
        File output_zarr  = "result/~{sample_id}_demux.zarr.zip"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/demuxlet:~{demuxlet_version}"
        zones: zones
        memory: "~{memory}G"
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(extra_disk_space + size(input_bam,"GB") +  size(ref_genotypes,"GB"))  + " HDD"
        cpu: num_cpu
        preemptible: preemptible
    }
}
