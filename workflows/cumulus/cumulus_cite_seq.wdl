version 1.0

workflow cumulus_cite_seq {
    input {
        File input_sample_sheet
        String output_directory
        String? genome
        # A CSV file containing the IgG control information for each antibody.
        File? antibody_control_csv

        String docker_registry = "cumulusprod"
        String cumulus_version = "0.15.0"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        Int memory = 10
        Int disk_space = 20
        Int preemptible = 2
    }

    call generate_cite_seq_tasks as Config {
        input:
            input_sample_sheet = input_sample_sheet,
            docker_registry = docker_registry,
            zones = zones,
            preemptible = preemptible
    }

    if (Config.cite_seq_ids[0] != '') {
        scatter (cite_seq_id in Config.cite_seq_ids) {
            call run_cumulus_merge_rna_adt {
                input:
                    output_name = cite_seq_id,
                    output_directory = output_directory,
                    input_raw_gene_bc_matrices_h5 = Config.id2rna[cite_seq_id],
                    input_adt_csv = Config.id2adt[cite_seq_id],
                    antibody_control_csv = antibody_control_csv,
                    docker_registry = docker_registry,
                    cumulus_version = cumulus_version,
                    zones = zones,
                    memory = memory,
                    disk_space = disk_space,
                    preemptible = preemptible
            }
        }
    }

    output {
        String output_folder = "~{output_directory}"
    }
}

task generate_cite_seq_tasks {
    input {
        File input_sample_sheet
        String docker_registry
        String zones
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import re
        import pandas as pd

        df = pd.read_csv('~{input_sample_sheet}', header = 0, dtype = str, index_col = False)
        for c in df.columns:
            df[c] = df[c].str.strip()
        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['OUTNAME'].str.contains(regex_pat)):
            print('OUTNAME must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)

        with open('cite_seq.txt', 'w') as fo_ids, open('id2rna.txt', 'w') as fo_rnas, open('id2adt.txt', 'w') as fo_adts:
            for idx, row in df.iterrows():
                fo_ids.write(row['OUTNAME'] + '\n')
                fo_rnas.write(row['OUTNAME'] + '\t' + row['RNA'] + '\n')
                fo_adts.write(row['OUTNAME'] + '\t' + row['ADT'] + '\n')

                if 'TYPE' in df.columns:
                    assert row['TYPE'] == 'cite-seq'
        CODE
    }

    output {
        Array[String] cite_seq_ids = read_lines('cite_seq.txt')
        Map[String, String] id2rna = read_map('id2rna.txt')
        Map[String, String] id2adt = read_map('id2adt.txt')
    }

    runtime {
        docker: "~{docker_registry}/config"
        zones: zones
        preemptible: "~{preemptible}"
    }
}

task run_cumulus_merge_rna_adt {
    input {
        String output_name
        String output_directory
        File input_raw_gene_bc_matrices_h5
        File input_adt_csv
        File? antibody_control_csv

        String docker_registry
        String cumulus_version
        String zones
        Int memory
        Int disk_space
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        from subprocess import check_call
        call_args = ['pegasus', 'merge_rna_adt', '~{input_raw_gene_bc_matrices_h5}', '~{input_adt_csv}', '~{output_name}']
        if '~{antibody_control_csv}' is not '':
            call_args.extend(['--antibody-control-csv', '~{antibody_control_csv}'])
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        gsutil -q cp ~{output_name}.h5sc ~{output_directory}/~{output_name}/
        # mkdir -p ~{output_directory}/~{output_name}
        # cp ~{output_name}.h5sc ~{output_directory}/~{output_name}/
	}

    output {
        String output_folder = "~{output_directory}/~{output_name}"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cumulus:~{cumulus_version}"
        zones: zones
        memory: "~{memory}G"
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
    }
}