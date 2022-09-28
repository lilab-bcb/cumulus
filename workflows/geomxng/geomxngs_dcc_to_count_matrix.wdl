version 1.0

workflow geomxngs_dcc_to_count_matrix {
    input {
        File ini
        File dcc_zip
        File dataset
        File acronym_file = "gs://regev-lab/resources/geomxngspipeline/pkc-alias.tsv"
        String pkc
        File lab_worksheet
        String docker_registry = "quay.io/cumulus"
        String docker_version = "1.0.0"
        String output_directory
        Int cpu = 1
        Int extra_disk_space = 5
        String memory = "8GB"
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
        String backend = "gcp"
        String aws_queue_arn = ""
    }
    Map[String, String] acronym2url = read_map(acronym_file)
    Boolean is_pkc_url = sub(pkc, "^.+\\.(zip|pkc)$", "URL") == "URL"
    File pkc_file = (if is_pkc_url then pkc else acronym2url[pkc])

    parameter_meta {
        dcc_zip:"DCC zip file from geomxngs_fastq_to_dcc workflow output"
        ini:"Configuration file in INI format, containing pipeline processing parameters"
        dataset:"Data QC and annotation file (Excel) downloaded from instrument after uploading DCC zip file; we only use the first tab (SegmentProperties)"
        pkc:"GeoMx DSP configuration file to associate assay targets with GeoMx HybCode barcodes and Seq Code primers"
        lab_worksheet:"A text file containing library setups"
        output_directory:"URL to write results (e.g. s3://foo/bar/out or gs://foo/bar/out)"
        backend:"Backend for computation (aws, gcp, or local)"
        docker_registry:"Docker registry"
        extra_disk_space:"Extra disk space in GB"
        preemptible: "Number of preemptible tries"
        memory:"Memory string"
        cpu:"Number of CPUs"
        aws_queue_arn:"The arn URI of the AWS job queue to be used (e.g. arn:aws:batch:us-east-1:xxxxx). Only works when backend is aws."
    }

    output {
        String count_matrix_h5ad = geomxngs_dcc_to_count_matrix_task.count_matrix_h5ad
        String count_matrix_text = geomxngs_dcc_to_count_matrix_task.count_matrix_text
        String count_matrix_metadata = geomxngs_dcc_to_count_matrix_task.count_matrix_metadata
    }

    call geomxngs_dcc_to_count_matrix_task {
        input:
            ini = ini,
            dcc_zip = dcc_zip,
            dataset = dataset,
            pkc = pkc_file,
            lab_worksheet = lab_worksheet,
            output_directory = output_directory,
            cpu = cpu,
            memory = memory,
            extra_disk_space = extra_disk_space,
            preemptible = preemptible,
            zones = zones,
            backend = backend,
            aws_queue_arn = aws_queue_arn,
            docker_registry = docker_registry,
            docker_version=docker_version
    }
}

task geomxngs_dcc_to_count_matrix_task {
    input {
        File ini
        File dcc_zip
        File dataset
        File pkc
        File lab_worksheet
        String output_directory
        String memory
        Int cpu
        Int preemptible
        Int extra_disk_space
        String zones
        String backend
        String aws_queue_arn
        String docker_registry
        String docker_version
    }
    Int disk_space = ceil(extra_disk_space + size(dcc_zip, 'GB')*2 + size(ini, 'GB') + size(dataset, 'GB') + size(pkc, 'GB') + size(lab_worksheet, 'GB'))
    String output_directory_trailing_slash = sub(output_directory, "[/\\s]+$", "") + '/'

    command {
        set -e
        monitor_script.sh > monitoring.log &

        mkdir -p dcc
        unzip -q -d dcc ~{dcc_zip}
        python /software/scripts/create-counts.py \
        --dcc dcc \
        --ini ~{ini} \
        --pkc ~{pkc} \
        --lab-worksheet ~{lab_worksheet} \
        --dataset ~{dataset} \
        --out results

        strato cp --backend ~{backend} -m -r results/* ~{output_directory_trailing_slash}
    }

    output {
        File monitoring_log = "monitoring.log"
        String count_matrix_h5ad = output_directory_trailing_slash + "counts.h5ad"
        String count_matrix_text = output_directory_trailing_slash + "counts.txt"
        String count_matrix_metadata = output_directory_trailing_slash + "metadata.txt"
    }

    runtime {
        docker: "~{docker_registry}/geomxngs_dcc_to_count_matrix:~{docker_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: cpu
        preemptible: preemptible
        queueArn: aws_queue_arn
    }

}
