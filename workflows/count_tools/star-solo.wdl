workflow starsolo {
    File input_sample_sheet
    File ref_gtf_file
    File ref_fasta_file
    String star_index_directory
    String output_directory
    Int? num_cpu = 32
    String? solo_type = 'CB_UMI_Simple'
    File solo_white_list
    Int solo_umi_length

    String? docker_registry = "cumulusprod"
    String? star_version = '2.7.3a'
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? memory = 32

    call star_genome_generate {
        input:
            input_gtf_file = ref_gtf_file,
            input_fasta_file = ref_fasta_file,
            output_directory = star_index_directory,
            num_cpu = num_cpu,
            docker_registry = docker_registry,
            star_version = star_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    call run_star_solo {
        input:
            input_sample_sheet = input_sample_sheet,
            solo_type = solo_type,
            white_list = solo_white_list,
            umi_length = solo_umi_length,
            num_cpu = num_cpu,
            genome_dir = star_genome_generate.star_genome_dir,
            output_directory = output_directory,
            docker_registry = docker_registry,
            star_version = star_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    output {
        File final_log = run_star_solo.final_log
        File summary_csv = run_star_solo.summary_csv
        File barcodes_stats = run_star_solo.barcodes_stats
        File genes_stats = run_star_solo.genes_stats
    }

}

task star_genome_generate {
    File input_gtf_file
    File input_fasta_file
    String output_directory
    Int num_cpu

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

        mkdir star
        STAR --runMode genomeGenerate --runThreadN ${num_cpu} --genomeDir star --genomeFastaFiles ${input_fasta_file} --sjdbGTFfile ${input_gtf_file} --genomeSAsparseD 3
        gsutil -m rsync -r star ${output_directory}
        # mkdir -p ${output_directory}
        # cp -r star/* ${output_directory}
    }

    output {
        File star_genome_dir = 'star'
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

task run_star_solo {
    File input_sample_sheet
    String solo_type
    File white_list
    Int umi_length
    Int num_cpu
    File genome_dir
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

        python <<CODE
        import os
        import pandas as pd
        from subprocess import check_call

        df = pd.read_csv("${input_sample_sheet}")
        r1_fastqs = df['R1_fastq'].values
        r2_fastqs = df['R2_fastq'].values

        call_args = ['STAR', '--soloType', '${solo_type}', '--soloCBwhitelist', '${white_list}', '--soloUMIlen', '${umi_length}', '--genomeDir', '${genome_dir}', '--runThreadN', '${num_cpu}']

        file_ext = os.path.splitext(r1_fastqs[0])
        if file_ext == '.gz':
            call_args.extend(['--readFilesCommand', 'zcat'])

        call_args.extend(['--readFilesIn', ','.join(r2_fastqs), ','.join(r1_fastqs)])

        print(' '.join(call_args))
        check_call(['mkdir', 'starsolo'])
        call_args.extend(['--outFileNamePrefix', 'starsolo/'])
        check_call(call_args)
        CODE

        gsutil -m rsync -r starsolo ${output_directory}
        # mkdir -p ${output_directory}
        # cp -r starsolo/* ${output_directory}

    }

    output {
        File final_log = 'starsolo/Log.final.out'
        File summary_csv = 'starsolo/Solo.out/Gene/Summary.csv'
        File barcodes_stats = 'starsolo/Solo.out/Barcodes.stats'
        File genes_stats = 'starsolo/Solo.out/Features.stats'
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