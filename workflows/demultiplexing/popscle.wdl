version 1.0

workflow popscle {
    input {
        # Sample ID
        String sample_id
        # Output directory (gs url + path)
        String output_directory
        # Link to a RNA count matrix in hdf5 format
        File input_rna
        # Link to a position-sorted BAM
        File input_bam
        # Link to reference genotype files
        File ref_genotypes
        # Only demultiplex cells/nuclei with at least <min_num_genes> expressed genes
        Int min_num_genes
        # Minimum mapping quality to consider (lower MQ will be ignored) [default: 20]
        Int? min_MQ
        # Minimum distance to the tail (lower will be ignored) [default: 0]
        Int? min_TD
        # Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB  [default: "CB"]
        String? tag_group
        # Tag representing UMIs. For 10x genomics, use UB  [default: "UB"]
        String? tag_UMI

        # demuxlet-specific input fields
        # FORMAT field to extract the genotype, likelihood, or posterior from
        String field = "GT"
        # Grid of alpha to search for [default: "0.1,0.2,0.3,0.4,0.5"]
        String? alpha

        # freemuxlet-specific input fields
        # Number of samples multiplexed together
        Int nsample
        # A comma-separated list of donor names for renaming clusters achieved by souporcell
        String? donor_rename

        # # Google cloud zones
        String zones
        # Number of CPUs to request for souporcell per pair
        Int num_cpu
        # Memory size (integer) in GB needed for souporcell per pair
        Int memory
        # Extra disk space (integer) in GB needed for popscle per pair
        Int extra_disk_space
        # Number of preemptible tries
        Int preemptible
        # Which docker registry to use
        String docker_registry
        # Popscle version to use
        String popscle_version
    }


    call popscle_task {
        input:
            sample_id = sample_id,
            output_directory = output_directory,
            input_rna = input_rna,
            input_bam = input_bam,
            ref_genotypes = ref_genotypes,
            donor_rename = donor_rename,
            min_num_genes = min_num_genes,
            min_MQ = min_MQ,
            alpha = alpha,
            min_TD = min_TD,
            tag_group = tag_group,
            tag_UMI = tag_UMI,
            field = field,
            nsample = nsample,
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
        Int min_num_genes
        Int nsample
        String field
        Int? min_MQ
        String? alpha
        Int? min_TD
        String? tag_group
        String? tag_UMI
        String? donor_rename

        String docker_registry
        String popscle_version
        Int num_cpu
        Int memory
        Int extra_disk_space
        Int preemptible
        String zones
    }

    Int disk_space = ceil(extra_disk_space + size(input_rna,"GB") + size(input_bam,"GB") + size(ref_genotypes, "GB"))
    String algorithm = if nsample == 0 then 'demuxlet' else 'freemuxlet'

    command {
        set -e
        monitor_script.sh > monitoring.log &

        mkdir result
        python /software/extract_barcodes_from_rna.py ~{input_rna} "~{sample_id}".barcodes.tsv ~{min_num_genes}

        python <<CODE
        from subprocess import check_call

        call_args = ['popscle', 'dsc-pileup', '--sam', '~{input_bam}', '--vcf', '~{ref_genotypes}', '--group-list', '~{sample_id}.barcodes.tsv', '--out', '~{sample_id}.plp']
        if '~{tag_group}' is not '':
            call_args.extend(['--tag-group', '~{tag_group}'])
        if '~{tag_UMI}' is not '':
            call_args.extend(['--tag-UMI', '~{tag_UMI}'])
        if '~{min_MQ}' is not '':
            call_args.extend(['--min-MQ', '~{min_MQ}'])
        if '~{min_TD}' is not '':
            call_args.extend(['--min-TD', '~{min_TD}'])

        print(' '.join(call_args))
        check_call(call_args)

        call_args = ['popscle', '~{algorithm}', '--plp', '~{sample_id}.plp', '--out', 'result/~{sample_id}']
        if '~{algorithm}' == 'demuxlet':
            call_args.extend(['--vcf', '~{ref_genotypes}'])
            if '~{field}' is not '':
                call_args.extend(['--field', '~{field}'])
            if '~{alpha}' is not '':
                alpha_list = '~{alpha}'.split(',')
                prefix_list = ['--alpha'] * len(alpha_list)
                alpha_args = list(sum(list(zip(prefix_list, alpha_list)), ()))
                call_args.extend(alpha_args)
        else:
            call_args.extend(['--nsample', '~{nsample}'])

        print(' '.join(call_args))
        check_call(call_args)

        cluster_result = 'result/~{sample_id}.best' if '~{algorithm}' == 'demuxlet' else 'result/~{sample_id}.clust1.samples.gz'
        call_args = ['python', '/software/generate_zarr.py', cluster_result, '~{input_rna}', 'result/~{sample_id}_demux.zarr.zip', '--ref-genotypes', '~{ref_genotypes}']
        if '~{algorithm}' == 'freemuxlet':
            call_args.extend(['--cluster-genotypes', 'result/~{sample_id}.clust1.vcf.gz'])
            if '~{donor_rename}' is not '':
                call_args.extend(['--donor-names', '~{donor_rename}'])

        print(' '.join(call_args))
        check_call(call_args)
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
