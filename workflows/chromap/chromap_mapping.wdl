version 1.0

workflow chromap_mapping {
    input {
        # Chromap version
        String chromap_version = "0.1.5"
        # Sample ID
        String sample_id
        # A comma-separated list of input FASTQs directories (urls)
        String input_fastqs_directories
        # Output directory, URL
        String output_directory

        # Keywords or a URL to a tar.gz file
        String genome
        # Index TSV file
        File acronym_file

        # R1 Fastq pattern, default as 10x atac/multiome format
        String read1_fastq_pattern = "_S*_L*_R1_001.fastq.gz"
        # R2 Fastq pattern, default as 10x atac/multiome format
        String read2_fastq_pattern = "_S*_L*_R3_001.fastq.gz"
        # Barcode Fastq pattern, default as 10x atac/multiome format
        String barcode_fastq_pattern = "_S*_L*_R2_001.fastq.gz"
        # Barcode whitelist
        File? barcode_whitelist
        # Read format
        String? read_format
        # Replace barcode
        File? barcode_translate

        # Preset option; available options: chip, hic, atac
        String preset = "atac"

        # Split alignment
        Boolean? split_alignment
        # Max edit distance
        Int? max_edit_dist_e
        # Min number of minimizers
        Int? min_num_minimizer_s
        # INT1[,INT2] Skip minimizers occuring > INT1 [500] times. INT2 [1000] is the threshold for a second round of seeding.
        String? ignore_minimizer_times_f
        # Max insert size, only for paired-end read mapping
        Int? max_insert_size_l
        # Min MAPQ in range [0, 60] for mappings to be output [30]
        Int? min_mapq_q
        # Skip mapping the reads of length less than Min read length
        Int? min_read_length
        # Trim adapters on 3’. This only works for paired-end reads.
        Boolean? trim_adaptors
        # Remove PCR duplicates
        Boolean? remove_pcr_duplicates
        # Remove PCR duplicates at bulk level for single cell data
        Boolean? remove_pcr_duplicates_at_bulk_level
        # Remove PCR duplicates at cell level for bulk data
        Boolean? remove_pcr_duplicates_at_cell_level
        # Perform Tn5 shift, only when --SAM is NOT set
        Boolean? tn5_shift
        # Low memory (use for big datasets)
        Boolean? low_mem
        # Max Hamming distance allowed to correct a barcode, max allowed 2
        Int? bc_error_threshold
        # Min probability to correct a barcode
        Float? bc_probability_threshold

        # Output mappings not in whitelist
        Boolean? output_mappings_not_in_whitelist
        # Output format; choices bed, tagalign, sam, pairs
        String? output_format

        # Customized chromsome order
        File? chr_order

        #Natural chromosome order for pairs flipping
        File? pairs_natural_chr_order

        # Number of cpus per chromap job
        Int num_cpu = 8
        # Memory string
        String memory = "64G"

        # Disk space in GB
        Int disk_space = 200
        # Which docker registry to use: quay.io/cumulus (default) or cumulusprod
        String docker_registry = "quay.io/cumulus"
        # Google cloud zones, default to "us-central1-b"
        String zones = "us-central1-b"
        # Backend
        String backend = "gcp"
        # Number of preemptible tries
        Int preemptible = 2
        # Arn string of AWS queue
        String awsQueueArn = ""
    }
    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")
    String docker_registry_stripped = sub(docker_registry, "/+$", "")

    Map[String, String] acronym2gsurl = read_map(acronym_file)
    Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"
    File genome_file = (if is_url then genome else acronym2gsurl[genome])

    call run_chromap {
        input:
            chromap_version = chromap_version,
            read1_fastq_pattern = read1_fastq_pattern,
            read2_fastq_pattern = read2_fastq_pattern,
            barcode_fastq_pattern = barcode_fastq_pattern,
            sample_id = sample_id,
            output_directory = output_directory_stripped,
            input_fastqs_directories = input_fastqs_directories,
            genome_file = genome_file,
            output_mappings_not_in_whitelist = output_mappings_not_in_whitelist,
            output_format = output_format,
            read_format = read_format,
            preset = preset,
            barcode_whitelist = barcode_whitelist,
            barcode_translate = barcode_translate,
            split_alignment = split_alignment,
            max_edit_dist_e = max_edit_dist_e,
            min_num_minimizer_s = min_num_minimizer_s,
            ignore_minimizer_times_f = ignore_minimizer_times_f,
            max_insert_size_l = max_insert_size_l,
            min_mapq_q = min_mapq_q,
            min_read_length = min_read_length,
            trim_adaptors = trim_adaptors,
            remove_pcr_duplicates = remove_pcr_duplicates,
            remove_pcr_duplicates_at_bulk_level = remove_pcr_duplicates_at_bulk_level,
            remove_pcr_duplicates_at_cell_level = remove_pcr_duplicates_at_cell_level,
            tn5_shift = tn5_shift,
            low_mem = low_mem,
            bc_error_threshold = bc_error_threshold,
            bc_probability_threshold = bc_probability_threshold,
            chr_order = chr_order,
            pairs_natural_chr_order = pairs_natural_chr_order,
            disk_space = disk_space,
            docker_registry = docker_registry_stripped,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            preemptible = preemptible,
            awsQueueArn = awsQueueArn,
            backend = backend
    }


    output {
        String output_aln_directory = run_chromap.output_aln_directory
        File monitoringLog = run_chromap.monitoringLog
    }

}

task run_chromap {
    input {
        String chromap_version
        String read1_fastq_pattern
        String read2_fastq_pattern
        String barcode_fastq_pattern
        String sample_id
        String output_directory
        String input_fastqs_directories
        File genome_file
        String? preset
        File? barcode_whitelist
        File? barcode_translate
        Boolean? output_mappings_not_in_whitelist
        String? output_format
        String? read_format
        Boolean? split_alignment
        Int? max_edit_dist_e
        Int? min_num_minimizer_s
        String? ignore_minimizer_times_f
        Int? max_insert_size_l
        Int? min_mapq_q
        Int? min_read_length
        Boolean? trim_adaptors
        Boolean? remove_pcr_duplicates
        Boolean? remove_pcr_duplicates_at_bulk_level
        Boolean? remove_pcr_duplicates_at_cell_level
        Boolean? tn5_shift
        Boolean? low_mem
        Int? bc_error_threshold
        Float? bc_probability_threshold
        File? chr_order
        File? pairs_natural_chr_order
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        String awsQueueArn
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &
        mkdir -p genome_dir
        tar xf ~{genome_file} -C genome_dir --strip-components 1

        python <<CODE
        import re
        import os
        from fnmatch import fnmatch
        from subprocess import check_call, CalledProcessError, DEVNULL, STDOUT
        import sys

        def set_up_input_fastq_files(l, folder, pattern):
            file_list = [f for f in os.listdir(folder) if fnmatch(f, "*" + pattern)]
            file_list.sort()
            for f in file_list:
                l.append(folder + '/' + f)

        r1_list = list()
        r2_list = list()
        if '~{preset}' == 'atac':
            index_list = list()

        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = re.sub('/+$', '', directory) # remove trailing slashes
            target = "~{sample_id}_" + str(i)
            try:
                call_args = ['strato', 'exists', '--backend', '~{backend}', directory + '/~{sample_id}/']
                print(' '.join(call_args))
                check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                call_args = ['strato', 'sync', '--backend', '~{backend}', '-m', directory + '/~{sample_id}', target]
                print(' '.join(call_args))
                check_call(call_args)

                set_up_input_fastq_files(r1_list, target, '~{read1_fastq_pattern}')
                set_up_input_fastq_files(r2_list, target, '~{read2_fastq_pattern}')

                if '~{preset}' == 'atac':
                    set_up_input_fastq_files(index_list, target, '~{barcode_fastq_pattern}')

            except CalledProcessError:
                if not os.path.exists(target):
                    os.mkdir(target)
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}~{read1_fastq_pattern}' , target + '/']
                print(' '.join(call_args))
                check_call(call_args)
                set_up_input_fastq_files(r1_list, target, '~{read1_fastq_pattern}')

                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}~{read2_fastq_pattern}' , target + '/']
                print(' '.join(call_args))
                check_call(call_args)
                set_up_input_fastq_files(r2_list, target, '~{read2_fastq_pattern}')

                if '~{preset}' == 'atac':
                    call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}~{barcode_fastq_pattern}' , target + '/']
                    print(' '.join(call_args))
                    check_call(call_args)
                    set_up_input_fastq_files(index_list, target, '~{barcode_fastq_pattern}')

        if '~{preset}' == 'atac':
            assert len(r1_list) == len(r2_list) ==  len(index_list)
            index_fq = ",".join(index_list)
        else:
            assert len(r1_list) == len(r2_list)

        read1_fq = ",".join(r1_list)
        read2_fq = ",".join(r2_list)

        call_args = ['chromap', '--preset', '~{preset}', '-r', 'genome_dir/ref.fa',
                     '-x', 'genome_dir/ref.index', '-1', read1_fq,
                     '-2', read2_fq]

        if '~{preset}' not in ['atac', 'hic', 'chip']:
            print('Choose from following preset options only: atac, chip or hic.')
            sys.exit(1)


        def check_gz_file(input_file):
            ### If file is gzipped, ungzip it
            root, ext = os.path.splitext(input_file)
            if  ext == '.gz':
                gunzip_command = ['gunzip', '-f', input_file]
                print(' '.join(gunzip_command))
                check_call(gunzip_command)
                return root
            return input_file


        if '~{preset}' == 'atac':
            if index_fq != '':
                call_args.extend(['-b', index_fq])
            if '~{barcode_translate}' != '':
                call_args.extend(['--barcode-translate', check_gz_file('~{barcode_translate}')])
            if '~{barcode_whitelist}' != '':
                call_args.extend(['--barcode-whitelist', check_gz_file('~{barcode_whitelist}')])
            out_file_suffix = '.bed'
        if '~{preset}' == 'chip':
            out_file_suffix = '.bed'
        if '~{preset}' == 'hic':
            out_file_suffix = '.pairs'

        if '~{output_format}' in ['bed','tagalign','sam','pairs'] and '~{output_format}' != '':
            if '~{output_format}' == 'tagalign':
                call_args.append('--TagAlign')
            if '~{output_format}' == 'bed':
                call_args.append('--BED')
            if '~{output_format}' == 'sam':
                call_args.append('--SAM')
            if '~{output_format}' == 'pairs':
                call_args.append('--pairs')
            out_file_suffix = '.' + '~{output_format}'

        out_file = '~{sample_id}' + out_file_suffix
        call_args.extend(['-o',out_file])

        if '~{output_mappings_not_in_whitelist}':
            call_args.append('--output-mappings-not-in-whitelist')

        if '~{split_alignment}':
            call_args.extend(['--split-alignment'])

        if '~{max_edit_dist_e}' != '':
            call_args.extend(['-e', '~{max_edit_dist_e}'])
        if '~{min_num_minimizer_s}' != '':
            call_args.extend(['-s', '~{min_num_minimizer_s}'])
        if '~{ignore_minimizer_times_f}' != '':
            call_args.extend(['-f', '~{ignore_minimizer_times_f}'])
        if '~{max_insert_size_l}' != '':
            call_args.extend(['-l', '~{max_insert_size_l}'])
        if '~{min_mapq_q}' != '':
            call_args.extend(['-q', '~{min_mapq_q}'])
        if '~{min_read_length}' != '':
            call_args.extend(['--min-read-length', '~{min_read_length}'])
        if '~{trim_adaptors}':
            call_args.append('--trim-adapters')
        if '~{remove_pcr_duplicates}':
            call_args.append('--remove-pcr-duplicates')
        if '~{remove_pcr_duplicates_at_bulk_level}':
            call_args.append('--remove-pcr-duplicates-at-bulk-level')
        if '~{remove_pcr_duplicates_at_cell_level}':
            call_args.append('--remove-pcr-duplicates-at-cell-level')
        if '~{read_format}' != '':
            call_args.extend(['--read-format', '~{read_format}'])
        if '~{tn5_shift}':
            call_args.append('--Tn5-shift')
        if '~{low_mem}':
            call_args.append('--low-mem')
        if '~{bc_error_threshold}' != '':
            call_args.extend(['--bc-error-threshold', '~{bc_error_threshold}'])
        if '~{bc_probability_threshold}' != '':
            call_args.extend(['--bc-probability-threshold', '~{bc_probability_threshold}'])
        if '~{num_cpu}' != '':
            call_args.extend(['-t', '~{num_cpu}'])
        if '~{chr_order}' != '':
            call_args.extend(['--chr-order', '~{chr_order}'])
        if '~{pairs_natural_chr_order}' != '':
            call_args.extend(['--pairs-natural-chr-order', '~{pairs_natural_chr_order}'])

        print(' '.join(call_args))
        check_call(call_args)

        call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', out_file, '~{output_directory}/~{sample_id}/']
        print(' '.join(call_args))
        check_call(call_args)

        CODE
    }

    output {
        String output_aln_directory = "~{output_directory}/~{sample_id}"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/chromap:~{chromap_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}
