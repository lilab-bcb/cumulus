workflow starsolo {
    File input_sample_sheet
    String genome_url
    String chemistry
    String output_directory
    Int num_cpu
    String star_version
    String? docker_registry = "cumulusprod"
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? memory = 32

    ## Determine solo_type
    String solo_type = if (chemistry == 'tenXV2' || chemistry == 'tenXV3') then 'CB_UMI_Simple' else 'Droplet'

    File wl_index_file = "gs://regev-lab/resources/count_tools/whitelist_index.tsv"
    # File wl_index_file = "whitelist_index.tsv"
    Map[String, String] wl_index2gsurl = read_map(wl_index_file)
    String whitelist_url = wl_index2gsurl[chemistry]


    call run_star_solo {
        input:
            input_sample_sheet = input_sample_sheet,
            solo_type = solo_type,
            chemistry = chemistry,
            num_cpu = num_cpu,
            star_version = star_version,
            genome_url = genome_url + '/starsolo.tar.gz',
            whitelist_url = whitelist_url,
            output_directory = output_directory,
            docker_registry = docker_registry,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    output {
        File monitoringLog = run_star_solo.monitoringLog
        String output_folder = run_star_solo.output_folder
    }

}

task run_star_solo {
    File input_sample_sheet
    String solo_type
    String chemistry
    Int num_cpu
    String genome_url
    String whitelist_url
    String output_directory
    String docker_registry
    String star_version
    Int disk_space
    String zones
    Int memory
    Int preemptible

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        gsutil -q -m cp ${genome_url} genome.tar.gz
        tar -zxvf genome.tar.gz
        rm genome.tar.gz

        gsutil -q -m cp ${whitelist_url} whitelist.txt.gz
        gunzip whitelist.txt.gz

        mkdir result
        
        python <<CODE
        import os
        import pandas as pd
        from subprocess import check_call

        df = pd.read_csv("${input_sample_sheet}")
        r1_fastqs = df['R1_fastq'].values
        r2_fastqs = df['R2_fastq'].values

        call_args = ['STAR', '--soloType', '${solo_type}', '--soloCBwhitelist', 'whitelist.txt', '--genomeDir', 'starsolo', '--runThreadN', '${num_cpu}']
        if '${chemistry}' is 'tenXV3':
            call_args.extend(['--soloUMIlen', '12'])

        file_ext = os.path.splitext(r1_fastqs[0])[-1]
        if file_ext == '.gz':
            call_args.extend(['--readFilesCommand', 'zcat'])

        call_args.extend(['--readFilesIn', ','.join(r2_fastqs), ','.join(r1_fastqs)])
        call_args.extend(['--outFileNamePrefix', 'result/'])

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        gsutil -q -m rsync -r result ${output_directory}
        # mkdir -p ${output_directory}
        # cp -r result/* ${output_directory}

    }

    output {
        File monitoringLog = 'monitoring.log'
        String output_folder = '${output_directory}'
    }

    runtime {
        docker: "${docker_registry}/starsolo:${star_version}"
        zones: zones
        memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${preemptible}"
    }
}