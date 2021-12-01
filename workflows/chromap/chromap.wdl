version 1.0

workflow chromap_mapping {
    input {
        # Chromap version
        String chromap_version = "0.1.3"
        # Sample ID
        String sample_id
        # A comma-separated list of input FASTQs directories (gs urls)
        String input_fastqs_directories
        # Output directory, gs URL
        String output_directory

        # Keywords or a URL to a tar.gz file
        String genome
        # Index TSV file 
        File acronym_file
        
        # Read1
        String read1 =  "R1"
        # Read2
        String read2 =  "R2"
        # Barcode file
        String barcode_idx = "I1"

        # Preset option
        String preset = "atac"
        # Barcode whitelist
        File? barcode_whitelist
        # Max edit distance
        Int? max_edit_dist
        # Min number of minimizers
        Int? num_minimizer
        # Ignore minimizers occuring these many times
        Int? ignore_minimizer_times 
        # Max insert size, only for paired-end read mapping
        Int? max_insert
        # Min MAPQ in range [0, 60] for mappings to be output [30]
        Int? mapq
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
        Boolean low_mem = true
        # Max Hamming distance allowed to correct a barcode, max allowed 2 
        Int? bc_error_threshold
        # Min probability to correct a barcode 
        Float? bc_probability_threshold
        # Num of threads for mapping
        Int? threads

        # Customized chromsome order
        File? chr_order
        # Output format
        String output_format = "BED"
        
        #Natural chromosome order for pairs flipping
        File? pairs_natural_chr_order

        # Number of cpus per chromap job
        Int num_cpu = 1
        # Memory string, e.g. 57.6G
        String memory = "50G"

        # Disk space in GB
        Int disk_space = 500
        # Which docker registry to use: quay.io/cumulus (default) or cumulusprod
        String docker_registry = "quay.io/cumulus"
        # Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Backend
        String backend = "gcp"
        # Number of preemptible tries
        Int preemptible = 2
        # Max number of retries for AWS instance
        Int awsMaxRetries = 5  
    }
    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")
    String docker_registry_stripped = sub(docker_registry, "/+$", "")

    Map[String, String] acronym2gsurl = read_map(acronym_file)
    Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"
    File genome_file = (if is_url then genome else acronym2gsurl[genome])

    call chromap {
        input:
            chromap_version = chromap_version,
            read1 = read1,
            read2 = read2,
            barcode_idx = barcode_idx,
            sample_id = sample_id,
            output_directory = output_directory_stripped,
            input_fastqs_directories = input_fastqs_directories,
            genome_file = genome_file,
            output_format = output_format,
            preset = preset,
            barcode_whitelist = barcode_whitelist,
            max_edit_dist = max_edit_dist,
            num_minimizer = num_minimizer,
            ignore_minimizer_times = ignore_minimizer_times,
            max_insert = max_insert,
            mapq = mapq,
            min_read_length = min_read_length,
            trim_adaptors = trim_adaptors,
            remove_pcr_duplicates = remove_pcr_duplicates,
            remove_pcr_duplicates_at_bulk_level = remove_pcr_duplicates_at_bulk_level,
            remove_pcr_duplicates_at_cell_level = remove_pcr_duplicates_at_cell_level,
            tn5_shift = tn5_shift,
            low_mem = low_mem,
            bc_error_threshold = bc_error_threshold,
            bc_probability_threshold = bc_probability_threshold,
            threads = threads,
            chr_order = chr_order,
            pairs_natural_chr_order = pairs_natural_chr_order,
            disk_space = disk_space,
            docker_registry = docker_registry_stripped,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend                                    
    }
        output {
            String output_aln_directory = chromap.output_aln_directory
            File monitoringLog = chromap.monitoringLog
        }

}

task chromap {
    input {
            String chromap_version
            String read1
            String read2
            String barcode_idx
            String sample_id
            String output_directory
            String input_fastqs_directories
            File genome_file
            String preset
            File? barcode_whitelist
            String? output_format
            Int? max_edit_dist
            Int? num_minimizer
            Int? ignore_minimizer_times
            Int? max_insert
            Int? mapq
            Int? min_read_length
            Boolean? trim_adaptors
            Boolean? remove_pcr_duplicates
            Boolean? remove_pcr_duplicates_at_bulk_level
            Boolean? remove_pcr_duplicates_at_cell_level
            Boolean? tn5_shift
            Boolean? low_mem
            Int? bc_error_threshold
            Float? bc_probability_threshold
            Int? threads
            File? chr_order
            File? pairs_natural_chr_order
            String docker_registry
            String zones
            Int num_cpu
            String memory
            Int disk_space
            Int preemptible
            Int awsMaxRetries
            String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &
        mkdir -p genome_dir
        tar xf ~{genome_file} -C genome_dir --strip-components 1   

        python <<CODE
        from re import match, sub
        import os
        from subprocess import check_call, CalledProcessError, DEVNULL, STDOUT
        from packaging import version
        
        fastqs = []
        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = sub('/+$', '', directory) # remove trailing slashes
            target = '~{sample_id}_' + str(i)
            try:
                call_args = ['strato', 'exists', '--backend', '~{backend}', directory + '/~{sample_id}/']
                print(' '.join(call_args))
                check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                call_args = ['strato', 'cp', '--backend','~{backend}','-r', '-m', directory + '/~{sample_id}', target]
                print(' '.join(call_args))
                check_call(call_args)
            except CalledProcessError:
                if not os.path.exists(target):
                    os.mkdir(target)
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}' + '_S*_L*_*_001.fastq.gz' , target]
                print(' '.join(call_args))
                check_call(call_args)
            fl=os.listdir(target)      
            fastqs.extend(fl)
        
        read1_fq = ",".join(list(filter(lambda k: '_~{read1}_' in k, fastqs)))
        read2_fq = ",".join(list(filter(lambda k: '_~{read2}_' in k, fastqs)))
        index_fq = ",".join(list(filter(lambda k: '_~{barcode_idx}_' in k, fastqs)))

        call_args = ['chromap', '--preset', '~{preset}', '-r', 'genome_dir/ref.fa', 
                     '-x', 'genome_dir/ref.index', '-1', read1_fq, 
                     '-2', read2_fq]

        if '~{preset}' == 'atac':
            call_args.extend['-b', index_fq]
            if '~{barcode_whitelist}' != '':
                call_args.extend['--barcode-whitelist', '~{barcode_whitelist}']

        if '~{output_format}' not in ['BED','BEDPE','TagAlign']:
            print('Choose output formats from BED, BEDPE or TagAlign. User chosen format ' +  '~{output_format}' + ' not available.' , file = sys.stderr)
            sys.exit(1)
        else:
            if '~{output_format}' == 'TagAlign':
                out_file = 'aln.tagAlign'
                call_args.extend(['--TagAlign', '-o', out_file])
            elif '~{output_format}' == 'BEDPE':
                out_file = 'aln.bedpe'
                call_args.extend(["--BEDPE", '-o', out_file])
            else:
                out_file = "aln.bed"
                call_args.extend(["--BED", '-o', out_file]) 

        if '~{max_edit_dist}' != '':
            call_args.extend(['-e', '~{max_edit_dist}'])
        if '~{num_minimizer}' != '':
            call_args.extend(['-s', '~{num_minimizer}'])
        if '~{ignore_minimizer_times}' != '':
            call_args.extend(['-f', '~{ignore_minimizer_times}'])
        if '~{max_insert}' != '':
            call_args.extend(['-l', '~{max_insert}'])
        if '~{mapq}' != '':
            call_args.extend(['-q', '~{mapq}'])
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
        if '~{tn5_shift}':
            call_args.append('--Tn5-shift')
        if '~{low_mem}':
            call_args.append('--low-mem')
        if '~{bc_error_threshold}' != '':
            call_args.extend(['--bc-error-threshold', '~{bc_error_threshold}'])
        if '~{bc_probability_threshold}' != '':
            call_args.extend(['--bc-probability-threshold', '~{bc_probability_threshold}'])
        if '~{threads}' != '':
            call_args.extend(['--threads', '~{threads}'])
        if '~{chr_order}' != '':
            call_args.extend(['--chr-order', '~{chr_order}'])
        if '~{pairs_natural_chr_order}' != '':
            call_args.extend(['--pairs-natural-chr-order', '~{pairs_natural_chr_order}'])

        print(' '.join(call_args))
        check_call(call_args)

        call_args = ['strato', 'sync', '--backend', '~{backend}', '-m', out_file, '~{output_directory}/~{sample_id}']
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
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
