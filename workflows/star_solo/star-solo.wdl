version 1.0

workflow starsolo {
    input {
        # Input TSV sample sheet describing metadata of each sample
        File input_tsv_file
        # Genome reference
        String genome
        # Chemistry, choosing from tenX_v3 (for 10X V3 chemistry), tenX_v2 (for 10X V2 chemistry), DropSeq, SeqWell and custom
        String chemistry
        # Cell barcode start position (1-based coordinate)
        Int? CBstart
        # Cell barcode length
        Int? CBlen
        # UMI start position (1-based coordinate)
        Int? UMIstart
        # UMI length
        Int? UMIlen
        # Cell barcode white list
        File? CBwhitelist
        # URL of output directory
        String output_directory
        # Number of CPUs to request per sample
        Int num_cpu = 32
        # STAR version to use. Currently only support 2.7.6a
        String star_version = "2.7.6a"
        # Docker registry, default to quay.io/cumulus 
        String docker_registry = "quay.io/cumulus"
        # Version of docker image to run configuration on the sample sheet
        String config_version = "0.2"
        # Disk space in GB needed per sample
        Int disk_space = 500
        # Number of maximum preemptible tries allowed
        Int preemptible = 2
        # Google cloud zones to consider for execution
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Memory size in GB needed for count per sample
        Int memory = 120
        # backend choose from "gcp", "aws", "local"
        String backend = "gcp"
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    File ref_index_file = if backend == "gcp" then "gs://regev-lab/resources/count_tools/ref_index.tsv" else (if backend == "aws" then "AWS URL" else "ref_index.tsv")
    Map[String, String] ref_index2gsurl = read_map(ref_index_file)
    String genome_url = if sub(genome, "^gs://.*", "gs://")=="gs://" then genome else ref_index2gsurl[genome] + '/starsolo.tar.gz' # Need to modify for AWS etc.

    File wl_index_file = if backend == "gcp" then "gs://regev-lab/resources/count_tools/whitelist_index.tsv" else (if backend == "aws" then "AWS URL" else "whitelist_index.tsv")
    Map[String, String] wl_index2gsurl = read_map(wl_index_file)
    String whitelist_url = wl_index2gsurl[chemistry]

    if (whitelist_url != 'null') {
        File whitelist = whitelist_url
    }


    call generate_count_config {
        input:
            input_tsv_file = input_tsv_file,
            docker_registry = docker_registry,
            zones = zones,
            version = config_version,
            preemptible = preemptible,
            backend = backend
    }

    if (generate_count_config.sample_ids[0] != '') {
        scatter (sample_id in generate_count_config.sample_ids) {
            call run_star_solo {
                input:
                    sample_id = sample_id,
                    r1_fastqs = generate_count_config.id2r1[sample_id],
                    r2_fastqs = generate_count_config.id2r2[sample_id],
                    chemistry = chemistry,
                    num_cpu = num_cpu,
                    star_version = star_version,
                    genome = genome_url,
                    CBstart = CBstart,
                    CBlen = CBlen,
                    UMIstart = UMIstart,
                    UMIlen = UMIlen,
                    CBwhitelist = CBwhitelist,
                    whitelist = whitelist,
                    output_directory = output_directory_stripped,
                    docker_registry = docker_registry,
                    disk_space = disk_space,
                    zones = zones,
                    memory = memory,
                    preemptible = preemptible,
                    backend = backend
            }
        }
    }


    output {
        String output_folder = output_directory
    }

}

task generate_count_config {
    input {
        File input_tsv_file
        String docker_registry
        String zones
        String version
        Int preemptible
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import re
        import pandas as pd
        from subprocess import check_call

        df = pd.read_csv('~{input_tsv_file}', sep = '\t', header = 0, dtype = str, index_col = False)
        for c in df.columns:
            df[c] = df[c].str.strip()

        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['Sample'].str.contains(regex_pat)):
            print('Sample must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)

        with open('sample_ids.txt', 'w') as fo1, open('sample_r1.tsv', 'w') as fo2, open('sample_r2.tsv', 'w') as fo3:
            for idx, row in df.iterrows():
                fo1.write(row['Sample'] + '\n')

                if 'Flowcells' in df.columns: # Fetch R1 and R2 fastqs automatically.
                    input_dir_list = list(map(lambda x: x.strip(), row['Flowcells'].split(',')))
                    r1_list = []
                    r2_list = []
                    for directory in input_dir_list:
                        directory = re.sub('/+$', '', directory)

                        call_args = ['gsutil', 'ls', directory]
                        # call_args = ['ls', directory]
                        with open('list_dir.txt', 'w') as tmp_fo:
                            check_call(call_args, stdout=tmp_fo)

                        with open('list_dir.txt', 'r') as tmp_fin:
                            f_list = tmp_fin.readlines()
                            f_list = list(map(lambda s: s.strip(), f_list))

                        r1_files = [f for f in f_list if re.match('.*_R1_.*.fastq.gz', f)]
                        r2_files = [f for f in f_list if re.match('.*_R2_.*.fastq.gz', f)]
                        r1_files.sort()
                        r2_files.sort()
                        # r1_files = list(map(lambda s: directory+'/'+s, r1_files))
                        # r2_files = list(map(lambda s: directory+'/'+s, r2_files))

                        r1_list.extend(r1_files)
                        r2_list.extend(r2_files)

                else:  # R1 and R2 fastqs specified in sample sheet.
                    r1_list = list(map(lambda s: s.strip(), row['R1'].split(',')))
                    r2_list = list(map(lambda s: s.strip(), row['R2'].split(',')))

                fo2.write(row['Sample'] + '\t' + ','.join(r1_list) + '\n')
                fo3.write(row['Sample'] + '\t' + ','.join(r2_list) + '\n')
        CODE
    }

    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Map[String, String] id2r1 = read_map('sample_r1.tsv')
        Map[String, String] id2r2 = read_map('sample_r2.tsv')
    }

    runtime {
        docker: "~{docker_registry}/config:~{version}"
        zones: zones
        preemptible: "~{preemptible}"
    }
}

task run_star_solo {
    input {
        String sample_id
        String r1_fastqs
        String r2_fastqs
        String chemistry
        Int num_cpu
        File genome
        Int? CBstart
        Int? CBlen
        Int? UMIstart
        Int? UMIlen
        File? CBwhitelist
        File? whitelist
        String output_directory
        String docker_registry
        String star_version
        Int disk_space
        String zones
        Int memory
        Int preemptible
        String backend
    }


    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        mkdir genome_ref
        tar -zxf "~{genome}" -C genome_ref --strip-components 1
        rm "~{genome}"

        mkdir result

        python <<CODE
        import os
        from subprocess import check_call

        r1_list = '~{r1_fastqs}'.split(',')
        r2_list = '~{r2_fastqs}'.split(',')
        assert len(r1_list) == len(r2_list)

        def rename_fastq(in_file, out_file):
            call_args = ['gsutil', '-q', '-m', 'cp', in_file, out_file]
            # call_args = ['cp', in_file, out_file]
            print(' '.join(call_args))
            check_call(call_args)

        for i in range(len(r1_list)):
            file_ext = '.fastq.gz' if os.path.splitext(r1_list[i])[-1] == '.gz' else '.fastq'
            rename_fastq(r1_list[i], 'R1_'+str(i)+file_ext)
            rename_fastq(r2_list[i], 'R2_'+str(i)+file_ext)

        call_args = ['STAR', '--soloType', 'CB_UMI_Simple', '--genomeDir', 'genome_ref', '--runThreadN', '~{num_cpu}', '--outSAMtype', 'BAM', 'Unsorted', '--outSAMheaderHD', '\\@HD', 'VN:1.4', 'SO:unsorted']

        white_list = 'None'
        if '~{whitelist}' != '':
            white_list = '~{whitelist}'
        else:
            if '~{chemistry}' == 'custom' and '~{CBwhitelist}' != '':
                white_list = '~{CBwhitelist}'

        if '~{chemistry}' == 'tenX_v3':
            call_args.extend(['--soloCBwhitelist', white_list, '--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '12'])
        elif '~{chemistry}' == 'tenX_v2':
            call_args.extend(['--soloCBwhitelist', white_list, '--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '10'])
        elif '~{chemistry}' in ['SeqWell', 'DropSeq']:
            call_args.extend(['--soloCBwhitelist', 'None', '--soloCBstart', '1', '--soloCBlen', '12', '--soloUMIstart', '13', '--soloUMIlen', '8'])
        elif '~{chemistry}' == 'custom':
            call_args.extend(['--soloCBwhitelist', white_list, '--soloCBstart', '~{CBstart}', '--soloCBlen', '~{CBlen}', '--soloUMIstart', '~{UMIstart}', '--soloUMIlen', '~{UMIlen}'])

        if file_ext == '.fastq.gz':
            call_args.extend(['--readFilesCommand', 'zcat'])

        def set_up_readfiles(prefix, n_files, f_ext):
            return ','.join([prefix+str(n)+f_ext for n in range(n_files)])

        call_args.extend(['--readFilesIn', set_up_readfiles('R2_', len(r2_list), file_ext), set_up_readfiles('R1_', len(r1_list), file_ext)])
        call_args.extend(['--outFileNamePrefix', 'result/~{sample_id}_'])

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        gsutil -q -m rsync -r result "~{output_directory}/~{sample_id}"
        # mkdir -p "~{output_directory}/~{sample_id}"
        # cp -r result/* "~{output_directory}/~{sample_id}"
    }

    output {
        File monitoringLog = 'monitoring.log'
        String output_folder = '~{output_directory}/~{sample_id}'
    }

    runtime {
        docker: "~{docker_registry}/starsolo:~{star_version}"
        zones: zones
        memory: "~{memory}G"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_cpu}"
        preemptible: "~{preemptible}"
    }
}
