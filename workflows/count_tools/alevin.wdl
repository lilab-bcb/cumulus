workflow alevin {
    String sample_id
    File r1_fastq
    File r2_fastq
    String genome_url
    String chemistry
    String output_directory
    Int? num_cpu = 32

    String? docker_registry = "cumulusprod"
    String? alevin_version = '1.1'
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? memory = 32

    File wl_index_file = "gs://regev-lab/resources/count_tools/whitelist_index.tsv"
    # File wl_index_file = "whitelist_index.tsv"
    Map[String, String] wl_index2gsurl = read_map(wl_index_file)
    String whitelist_url = wl_index2gsurl[chemistry]

    String library_type = if (chemistry == 'tenXV2' || chemistry == 'tenXV3' || chemistry == 'Dropseq') then 'ISR' else ''

    call run_alevin {
        input:
            sample_id = sample_id,
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            genome = genome_url,
            chemistry = chemistry,
            library_type = library_type,
            whitelist = whitelist_url,
            num_cpu = num_cpu,
            output_directory = output_directory,
            docker_registry = docker_registry,
            alevin_version = alevin_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    output {
        File monitoringLog = run_alevin.monitoringLog
        String output_folder = run_alevin.output_folder
    }

}


task run_alevin {
    String sample_id
    File r1_fastq
    File r2_fastq
    File genome
    String library_type
    String chemistry
    File whitelist
    Int num_cpu
    String output_directory

    String docker_registry
    String alevin_version
    Int disk_space
    String zones
    Int memory
    Int preemptible

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        tar -zxvf ${genome}
        rm ${genome}

        python <<CODE
        from subprocess import check_call

        call_args = ['salmon', 'alevin', '-l', '${library_type}', '-1']
        call_args.extend(['-1', '${r1_fastq}'])
        call_args.extend(['-2', '${r2_fastq}'])
        call_args.extend(['-i', 'alevin-ref/salmon_index'])
        
        if '${chemistry}' is 'tenX_v3':
            call_args.append('--chromiumV3')
        elif '${chemistry}' is 'tenX_v2':
            call_args.append('--chromium')
        else:
            call_args.append('--dropseq')

        if '${chemistry}' in ['tenX_v2', 'tenX_v3']:
            call_args.extend(['--whitelist', '${whitelist}'])

        call_args.extend(['-p', '${num_cpu}', '-o', 'result', '--tgMap', 'alevin-ref/txp2gene.tsv'])

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        gsutil -q -m rsync -r result ${output_directory}/${sample_id}
        # mkdir -p ${output_directory}/${sample_id}
        # cp -r alevin_output/* ${output_directory}/${sample_id}
    }

    output {
        File monitoringLog = 'monitoring.log'
        String output_folder = '${output_directory}/${sample_id}'
    }

    runtime {
        docker: "${docker_registry}/alevin:${alevin_version}"
        zones: zones
        memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${preemptible}"
    }
}