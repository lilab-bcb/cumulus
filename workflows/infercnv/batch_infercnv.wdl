version 1.0

import "https://raw.githubusercontent.com/lilab-bcb/cumulus/infercnv/workflows/infercnv/infercnv.wdl" as icv

workflow batch_infercnv {
    input {
        File input_csv_file
        String output_directory
        String gene_ordering

        String acronym_file = "gs://cumulus-ref/resources/infercnv/index.tsv"
        String awsQueueArn = ""
    }

    String zones = "us-west1-a us-west1-b us-west1-c"
    Int preemptible = 2

    call generate_config {
        input:
            input_csv_file = input_csv_file,
            zones = zones,
            preemptible = preemptible
    }

    if (length(generate_config.sample_ids) > 0) {
        scatter (sample_id in generate_config.sample_ids) {
            call icv.InferCNV as InferCNV {
                input:
                    input_zarr = generate_config.sample2zarr[sample_id],
                    sample_id = sample_id,
                    output_directory = output_directory,
                    gene_ordering = gene_ordering,
                    acronym_file = acronym_file,
                    ref_group_names = generate_config.sample2ref[sample_id],
                    zones = zones,
                    preemptible = preemptible,
                    awsQueueArn = awsQueueArn
            }
        }
    }

    output {
        Array[String]? output_folders = InferCNV.output_folder
    }
}

task generate_config {
    input {
        File input_csv_file
        String zones
        Int preemptible
    }

    String docker_registry = "quay.io/cumulus"
    String config_version = "0.3"

    command <<<
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import pandas as pd

        df = pd.read_csv("~{input_csv_file}", header=0, dtype=str, index_col=False)
        df.columns = df.columns.str.strip()
        for c in df.columns:
            df[c] = df[c].str.strip()

        with open("sample_ids.txt", 'w') as fout_id, open("sample2zarr.txt", 'w') as fout_zarr, open("sample2ref.txt", 'w') as fout_ref:
            for idx, row in df.iterrows():
                fout_id.write(row['Sample']+'\n')
                fout_zarr.write(row['Sample']+'\t'+row['Zarr']+'\n')
                fout_ref.write(row['Sample']+'\t'+row['RefGroups']+'\n')
        CODE
    >>>

    output {
        Array[String] sample_ids = read_lines("sample_ids.txt")
        Map[String, String] sample2zarr = read_map("sample2zarr.txt")
        Map[String, String] sample2ref = read_map("sample2ref.txt")
    }

    runtime {
        docker: "~{docker_registry}/config:~{config_version}"
        zones: zones
        preemptible: preemptible
    }
}
