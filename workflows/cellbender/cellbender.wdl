version 1.0

workflow cellbender {
  input {
      # Input csv file
      File input_csv_file
      # Outputs
      String output_directory
      # Docker image for cellbender remove-background version
      String docker_registry = "quay.io/cumulus"
      # Cellbender version to use. Currently support: 0.2.0
      String cellbender_version = "0.2.0"      

      # Expected cells
      Int? expected_cells
      # The number of droplets from the rank-ordered UMI plot that will be analyzed
      Int? total_droplets_included
      # Which model is being used for count data. ["simple", "ambient", "swapping", "full"]
      String? model
      # Droplets with UMI counts below this number are completely excluded from the analysis
      Int? low_count_threshold
      # Target false positive rate in (0, 1).  A false positive is a true signal count that is erroneously removed.
      String? fpr
      # Number of epochs to train.
      Int? epochs
      # Dimension of latent variable z.
      Int? z_dim
      # Dimension of hidden layers in the encoder for z.
      String? z_layers
      # Training detail: the fraction of the training data each epoch that is drawn (randomly sampled) from surely empty droplets.
      Float? empty_drop_training_fraction
      # Integer indices of genes to ignore entirely.
      String? blacklist_genes
      # Training detail: lower learning rate for inference.
      Float? learning_rate
      # Cause remove-background to operate on gene counts only, ignoring other features.
      Boolean exclude_antibody_capture = false

      # Google cloud zones to consider for execution
      String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
      # Disk space in GB needed per sample
      Int disk_space = 50
      # Boot disk space in GB needed per sample
      Int boot_disk_size_GB = 20
      # Number of maximum preemptible tries allowed
      Int preemptible = 0
      # Number of maximum retries when running on AWS
      Int awsMaxRetries = 5
      # Arn string of AWS queue
      String awsQueueArn = ""
      # Number of CPUs to request per sample
      Int num_cpu = 4
      # Memory size string represent memory needed for count per sample
      String memory = "15G"
      # GPU Type
      String gpu_type = "nvidia-tesla-t4"
      # backend choose from "gcp", "aws", "local"
      String backend = "gcp"
  }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

  call generate_config {
      input:
          input_csv_file = input_csv_file,
          docker_registry = docker_registry,
          zones = zones,
          cellbender_version = cellbender_version,
          preemptible = preemptible,
          awsMaxRetries = awsMaxRetries,
          backend = backend
    } 

  if (length(generate_config.sample_ids) > 0) {
        scatter (sample_id in generate_config.sample_ids) {
            call run_cellbender_remove_background_gpu {
                input:
                    sample_name = sample_id,
                    input_10x_h5_file = generate_config.sample2dir[sample_id],
                    output_directory = output_directory_stripped,
                    docker_registry = docker_registry,
                    expected_cells = expected_cells,
                    total_droplets_included = total_droplets_included,
                    model = model,
                    low_count_threshold = low_count_threshold,
                    fpr = fpr,
                    epochs = epochs,
                    z_dim = z_dim,
                    z_layers = z_layers,
                    empty_drop_training_fraction = empty_drop_training_fraction,
                    blacklist_genes = blacklist_genes,
                    learning_rate = learning_rate,
                    exclude_antibody_capture = exclude_antibody_capture,
                    zones = zones,
                    disk_space = disk_space,
                    cellbender_version = cellbender_version,
                    boot_disk_size_GB = boot_disk_size_GB,
                    num_cpu = num_cpu,
                    memory = memory,
                    preemptible = preemptible,
                    awsMaxRetries = awsMaxRetries,
                    awsQueueArn = awsQueueArn,
                    gpu_type = gpu_type,
                    backend = backend
            }
        }
    }        

  output {
    Array[String]? count_outputs = run_cellbender_remove_background_gpu.output_dir
    Array[File]? monitoringLog = run_cellbender_remove_background_gpu.monitoringLog
  }

}    

task generate_config {
    input {
        File input_csv_file
        String zones
        String cellbender_version
        String docker_registry
        Int preemptible
        Int awsMaxRetries
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import re, sys
        import pandas as pd
        from subprocess import check_call

        df = pd.read_csv('~{input_csv_file}', header = 0, dtype = str, index_col = False, sep="\t")
        df.columns = df.columns.str.strip()
        for c in df.columns:
            df[c] = df[c].str.strip()
        
        cellbender_cols = ["Sample","Location"]
        for c in cellbender_cols:
            if c not in df.columns:
                print("To run Cellbender, following column is required: " + c + ". Please correct sample sheet.", file = sys.stderr)
                sys.exit(1)       

        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['Sample'].str.contains(regex_pat)):
            print('Sample must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)

        with open('sample_ids.txt', 'w') as fo1, open('sample2dir.txt', 'w') as fo2:
            for sample_id in df['Sample'].unique():
                df_local = df.loc[df['Sample'] == sample_id]
                dirs = df_local['Location'].values
                fo1.write(sample_id + '\n')
                fo2.write(sample_id + '\t' + ','.join(dirs) + '\n')
        CODE
    }
    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Map[String, String] sample2dir = read_map('sample2dir.txt')
    }

    runtime {
        docker: "~{docker_registry}/cellbender:~{cellbender_version}"
        zones: zones
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }    
}

task run_cellbender_remove_background_gpu {
    input {
        String sample_name
        File input_10x_h5_file
        String output_directory
        String docker_registry
        String cellbender_version
        Int? expected_cells
        Int? total_droplets_included
        String? model
        Int? low_count_threshold
        String? fpr
        Int? epochs
        Int? z_dim
        String? z_layers
        Float? empty_drop_training_fraction
        String? blacklist_genes
        Float? learning_rate
        Boolean exclude_antibody_capture
        String zones
        Int disk_space
        Int boot_disk_size_GB
        Int preemptible
        Int awsMaxRetries
        String awsQueueArn
        Int num_cpu
        String memory
        String gpu_type
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        cellbender remove-background \
           --input "~{input_10x_h5_file}" \
           --output "~{sample_name}_out.h5" \
           --cuda \
           ~{"--expected-cells " + expected_cells} \
           ~{"--total-droplets-included " + total_droplets_included} \
           ~{"--fpr " + fpr} \
           ~{"--model " + model} \
           ~{"--low-count-threshold " + low_count_threshold} \
           ~{"--epochs " + epochs} \
           ~{"--z-dim " + z_dim} \
           ~{"--z-layers " + z_layers} \
           ~{"--empty-drop-training-fraction " + empty_drop_training_fraction} \
           ~{"--blacklist-genes " + blacklist_genes} \
           ~{"--learning-rate " + learning_rate} \
           ~{true="--exclude-antibody-capture" false=" " exclude_antibody_capture}

        strato cp --backend ~{backend} ~{sample_name}_out* ~{output_directory}/${sample_name}/
    }

    output {
        File monitoringLog = "monitoring.log"
        String output_dir = "~{output_directory}/${sample_name}"
    }

    runtime {
        docker: "~{docker_registry}/cellbender:~{cellbender_version}"
        bootDiskSizeGb: "~{boot_disk_size_GB}"
        disks: "local-disk ~{disk_space} HDD"
        memory: memory
        cpu: num_cpu
        zones: zones
        gpuCount: 1
        gpuType: "~{gpu_type}"
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
        queueArn: awsQueueArn
    }

}

