version 1.0

workflow starsolo {
    input {
        String sample_id
        File r1_fastq
        File r2_fastq
        String genome_url
        String chemistry
        String output_directory
        Int num_cpu = 32
        String star_version
        String docker_registry = "cumulusprod"
        Int disk_space = 500
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        Int memory = 120
    }

    ## Determine solo_type
    String solo_type = if (chemistry == 'tenX_v2' || chemistry == 'tenX_v3') then 'CB_UMI_Simple' else 'Droplet'

    File wl_index_file = "gs://regev-lab/resources/count_tools/whitelist_index.tsv"
    # File wl_index_file = "whitelist_index.tsv"
    Map[String, String] wl_index2gsurl = read_map(wl_index_file)
    String whitelist_url = wl_index2gsurl[chemistry]


    call run_star_solo {
        input:
            sample_id = sample_id,
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            solo_type = solo_type,
            chemistry = chemistry,
            num_cpu = num_cpu,
            star_version = star_version,
            genome = genome_url,
            whitelist = whitelist_url,
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
    input {
        String sample_id
        File r1_fastq
        File r2_fastq
        String solo_type
        String chemistry
        Int num_cpu
        File genome
        File whitelist
        String output_directory
        String docker_registry
        String star_version
        Int disk_space
        String zones
        Int memory
        Int preemptible
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

        call_args = ['STAR', '--soloType', '~{solo_type}', '--soloCBwhitelist', '~{whitelist}', '--genomeDir', 'genome_ref', '--runThreadN', '~{num_cpu}', '--outSAMtype', 'BAM', 'Unsorted', '--outSAMheaderHD', '\\@HD', 'VN:1.4', 'SO:unsorted']
        if '~{chemistry}' is 'tenX_v3':
            call_args.extend(['--soloUMIlen', '12'])

        file_ext = os.path.splitext('~{r1_fastq}')[-1]
        if file_ext == '.gz':
            call_args.extend(['--readFilesCommand', 'zcat'])

        call_args.extend(['--readFilesIn', '~{r2_fastq}', '~{r1_fastq}'])
        call_args.extend(['--outFileNamePrefix', 'result/~{sample_id}_'])

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        gsutil -q -m rsync -r result ~{output_directory}/~{sample_id}
        # mkdir -p ~{output_directory}/~{sample_id}
        # cp -r result/* ~{output_directory}/~{sample_id}

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