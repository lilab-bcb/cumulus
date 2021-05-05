version 1.0

workflow popscle {
    input {
        String sample_id
        String output_directory
        File input_rna
        File input_bam

        Int min_num_genes = 100
        String field = "GT"
        Int? min_MQ
        Int? min_TD


        File ref_genotypes

        # demuxlet
        String? alpha
        Float? geno_error

        # freemuxlet
        Int? nsample = 4

        String? tag_group = "CB"
        String? tag_UMI = "UB"
        String zones = "us-central1-b us-east1-d us-west1-a us-west1-b"
        Int num_cpu = 1
        Int memory = 10
        Int extra_disk_space = 2
        Int preemptible = 2
        String docker_registry = "quay.io/cumulus"
        String popscle_version = "0.1b"
        String algorithm
    }



    call popscle_task {
        input:
            sample_id = sample_id,
            output_directory = output_directory,
            input_rna = input_rna,
            input_bam = input_bam,
            ref_genotypes = ref_genotypes,
            algorithm = algorithm,
            min_num_genes = min_num_genes,
            min_MQ = min_MQ,
            alpha = alpha,
            min_TD = min_TD,
            tag_group = tag_group,
            tag_UMI = tag_UMI,
            field = field,
            nsample = nsample,
            geno_error = geno_error,
            docker_registry = docker_registry,
            popscle_version = popscle_version,
            num_cpu = num_cpu,
            memory = memory,
            extra_disk_space = extra_disk_space,
            zones = zones,
            preemptible = preemptible
    }

    output {
        String output_folder = popscle_task.output_folder
        File output_zarr = popscle_task.output_zarr
        File monitoringLog = popscle_task.monitoringLog
    }

}

task popscle_task {
    input {
        String sample_id
        String output_directory
        File input_rna
        File input_bam
        File ref_genotypes
        String algorithm
        Int min_num_genes
        Int? nsample
        String field
        Int? min_MQ
        String? alpha
        Int? min_TD
        String? tag_group
        String? tag_UMI
        Float? geno_error

        String docker_registry
        String popscle_version
        Int num_cpu
        Int memory
        Int extra_disk_space
        Int preemptible
        String zones
    }

    String type  = if algorithm=='demuxlet' then "demuxlet" else "dsc-pileup"
    Int disk_space = ceil(extra_disk_space + size(input_rna,"GB") + size(input_bam,"GB") + size(ref_genotypes, "GB"))

    command {
        set -e
        monitor_script.sh > monitoring.log &

        mkdir result
        python /software/extract_barcodes_from_rna.py ~{input_rna} "~{sample_id}".barcodes.tsv ~{min_num_genes}

        python <<CODE
        from subprocess import check_call
        type = '~{type}'


        call_args = ['popscle', '~{type}',  '--sam', '~{input_bam}', '--group-list', '~{sample_id}.barcodes.tsv', '--out', 'result/~{sample_id}', '--vcf', '~{ref_genotypes}']

        if type == 'demuxlet':
            call_args.extend(['--field', '~{field}'])

        if '~{min_MQ}' is not '':
            call_args.extend(['--min-MQ', '~{min_MQ}'])
        if '~{min_TD}' is not '':
            call_args.extend(['--min-TD', '~{min_TD}'])
        if '~{tag_group}' is not '':
            call_args.extend(['--tag-group', '~{tag_group}'])
        if '~{tag_UMI}' is not '':
            call_args.extend(['--tag-UMI', '~{tag_UMI}'])
        if type == 'demuxlet':
            if '~{alpha}' is not '':
                alpha_list = '~{alpha}'.split(',')
                prefix_list = ['--alpha'] * len(alpha_list)
                alpha_args = list(sum(list(zip(prefix_list, alpha_list)), ()))
                call_args.extend(alpha_args)
            if '~{geno_error}' is not '':
                call_args.extend(['--geno-error-offset', '~{geno_error}'])

        print(' '.join(call_args))
        check_call(call_args)

        if type != 'demuxlet':
            call_args = ['popscle', 'freemuxlet', '--plp', 'result/~{sample_id}', '--out', 'result/~{sample_id}']
            if '~{nsample}' is not '':
                call_args.extend(['--nsample', '~{nsample}'])
            print(' '.join(call_args))
            check_call(call_args)
            check_call(['python', '/software/generate_zarr.py', 'result/~{sample_id}.clust1.samples.gz', '~{input_rna}', 'result/~{sample_id}_demux.zarr.zip'])
        else:
            check_call(['python', '/software/generate_zarr.py', 'result/~{sample_id}.best', '~{input_rna}', 'result/~{sample_id}_demux.zarr.zip'])
        CODE

        gsutil -q -m rsync -r result "~{output_directory}/~{sample_id}"

        # mkdir -p "~{output_directory}/~{sample_id}"
        # cp result/* "~{output_directory}/~{sample_id}"
    }

    output {
        String output_folder = "~{output_directory}/~{sample_id}"
        File output_zarr  = "result/~{sample_id}_demux.zarr.zip"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/popscle:~{popscle_version}"
        zones: zones
        memory: "~{memory}G"
        bootDiskSizeGb: 12
        disks: "local-disk " + disk_space  + " HDD"
        cpu: num_cpu
        preemptible: preemptible
    }
}
