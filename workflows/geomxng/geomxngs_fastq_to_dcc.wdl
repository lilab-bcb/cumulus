version 1.0

workflow geomxngs_fastq_to_dcc {
    input {
        File ini
        String fastq_directory
        String output_directory
        File? fastq_rename
        String docker_registry = "quay.io/cumulus"
        String geomxngs_version = "2.3.3.10"
        Int cpu = 4
        Int disk_space = 500
        String memory = "64GB"
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
        String aws_queue_arn = ""
        Boolean delete_fastq_directory = false
    }

    parameter_meta {
        ini:"Configuration file in INI format, containing pipeline processing parameters"
        fastq_directory:"FASTQ directory URL (e.g. s3://foo/bar/fastqs or gs://foo/bar/fastqs). Separate multiple directories with a comma"
        output_directory:"URL to write results (e.g. s3://foo/bar/out or gs://foo/bar/out)"
        fastq_rename:"Optional 2 column TSV file with no header used to map original FASTQ names to FASTQ names that GeoMX recognizes"
        docker_registry :"Docker registry"
        geomxngs_version : "Version of the geomx software, currently only 2.3.3.10"
        cpu:"Number of CPUs"
        disk_space : "Disk space in GB"
        memory : "Memory string"
        preemptible :"Number of preemptible tries"
        zones : "Google cloud zones"
        aws_queue_arn : "The arn URI of the AWS job queue to be used (e.g. arn:aws:batch:us-east-1:xxxxx). Only works when backend is aws"
        delete_fastq_directory:"Whether to delete the input fastqs upon successful completion"
    }

    output {
        String geomxngs_output = geomxngs_task.geomxngs_output
        String dcc_zip = geomxngs_task.dcc_zip
        File local_dcc_zip = geomxngs_task.local_dcc_zip
    }

    call geomxngs_task {
        input:
            fastq_directory = fastq_directory,
            delete_fastq_directory=delete_fastq_directory,
            ini = ini,
            output_directory = output_directory,
            fastq_rename = fastq_rename,
            docker_registry = docker_registry,
            geomxngs_version = geomxngs_version,
            cpu = cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            zones = zones,
            aws_queue_arn = aws_queue_arn
    }

}

task geomxngs_task {
    input {
        String fastq_directory
        String output_directory
        File ini
        File? fastq_rename
        String memory
        Int cpu
        Int preemptible
        Int disk_space
        String docker_registry
        String geomxngs_version
        String zones
        String aws_queue_arn
        Boolean delete_fastq_directory
    }
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    command {
        set -e
        export DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1
        monitor_script.sh > monitoring.log &

        mkdir fastqs

        python <<CODE
        import os
        import pandas as pd
        from subprocess import check_call
        import configparser
        # update cpus in ini file
        ini_path = '~{ini}'
        cpu = '~{cpu}' # can only write strings to ini file
        output_path = 'local.ini'

        config = configparser.ConfigParser()
        config.optionxform = str  # prevent conversion of keys to lowercase
        config.read(ini_path)
        processing_keys = ["Processing", "Processing_v2"]
        found = False
        for processing_key in processing_keys:
            if processing_key in config:
                config[processing_key]["threads"] = cpu
                found = True
                break
        if not found:
            raise ValueError("Processing section not found")
        with open(output_path, "wt") as out:
            config.write(out)

        # download and rename fastqs
        rename = '~{fastq_rename}'
        remote_fastq_dirs = '~{fastq_directory}'.split(',')
        local_fastq_dirs = []
        if len(remote_fastq_dirs) == 1:
            local_fastq_dirs.append('fastqs')
            check_call(['strato', 'sync', '-m', remote_fastq_dirs[0], 'fastqs/'])
        else:  # geomx pipeline only works with one directory of fastqs
            for i in range(len(remote_fastq_dirs)):
                local_fastq_dir = 'fastqs-' + str(i + 1)
                local_fastq_dirs.append(local_fastq_dir)
                os.makedirs(local_fastq_dir, exist_ok=True)
                check_call(
                    ['strato', 'sync', '-m', remote_fastq_dirs[i], local_fastq_dir])

        if rename:
            df = pd.read_csv(rename, sep="\t", header=None, names=["original_name", "new_name"])
            df = df.dropna()
            # strip path
            df["original_name"] = df["original_name"].str.split("/").str[-1]
            df["new_name"] = df["new_name"].str.split("/").str[-1]

            for i in range(len(df)):
                d = df.iloc[i]
                original_name = d["original_name"]
                new_name = d["new_name"]
                for local_fastq_dir in local_fastq_dirs:
                    src = os.path.join(local_fastq_dir, original_name)
                    if os.path.exists(src):
                        os.rename(src, os.path.join(local_fastq_dir, new_name))
        # Illumina convention: SampleName_S1_L001_R1_001.fastq.gz
        # add a digit to sample name to ensure it's unique. For example:
        # HITS5483936 to HITS5483936-1
        if len(local_fastq_dirs) > 1:
            local_dir = "fastqs"
            for local_fastq_dir in local_fastq_dirs:
                for f in os.listdir(local_fastq_dir):
                    file_name = os.path.basename(f)
                    dest = os.path.join(local_dir, file_name)
                    counter = 1
                    while os.path.exists(dest):
                        name_tokens = file_name.split('_')
                        name_tokens[1] = name_tokens[1] + '-' + str(counter)
                        dest = os.path.join(local_dir, '_'.join(name_tokens))
                        counter = counter + 1
                    os.rename(os.path.join(local_fastq_dir, f), dest)
        CODE

        export TMPDIR="/tmp"
        geomx_expect.exp
        strato sync -m results ~{output_directory_stripped}

        python <<CODE
        from subprocess import check_call
        delete_fastq_directory = '~{delete_fastq_directory}' == 'true'
        if delete_fastq_directory:
            urls = '~{fastq_directory}'.split(',')
            for url in urls:
                if not url.endswith("/"):
                    url += "/"
                try:
                    call_args = ["strato", "rm", "-m", "-r", url]
                    check_call(call_args)
                    print("Deleted " + url)
                except:
                    print("Failed to delete " + url)
        CODE
    }

    output {
        File monitoring_log = "monitoring.log"
        String geomxngs_output = "${output_directory_stripped}"
        String local_dcc_zip = glob("results/*.zip")[0]
        String dcc_zip = output_directory_stripped  + "/" + basename(local_dcc_zip)
    }

    runtime {
        docker: "~{docker_registry}/geomxngs_fastq_to_dcc:~{geomxngs_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: cpu
        preemptible: preemptible
        queueArn: aws_queue_arn
    }

}
