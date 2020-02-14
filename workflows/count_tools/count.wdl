version 1.0

import "merge_fastqs.wdl" as mfs
import "star-solo.wdl" as sts
import "bustools.wdl" as kbc
import "alevin.wdl" as ale
import "optimus-count.wdl" as opm

workflow count {
    input {
        File input_tsv_file

        String genome
        String chemistry
        String output_directory

        # Count
        Boolean run_count = true 
        String count_tool

        String docker_registry = "cumulusprod"
        Int num_cpu = 32
        Int disk_space = 100
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        Int memory = 32

        # Star-Solo
        String starsolo_star_version = "2.7.3a"

        # Alevin
        String alevin_version = '1.1'

        # Bustools
        Boolean bustools_output_loom = false
        Boolean bustools_output_h5ad = false
        String bustools_docker = "shaleklab/kallisto-bustools"
        String bustools_version = '0.24.4'

        # Optimus
        String optimus_version = 'optimus_v1.4.0'
        Boolean optimus_output_loom = false
    }

    call generate_count_config {
        input:
            input_tsv_file = input_tsv_file,
            docker_registry = docker_registry,
            zones = zones,
            preemptible = preemptible
    }
    
    File ref_index_file = "gs://regev-lab/resources/count_tools/ref_index.tsv"
    # File ref_index_file = "ref_index.tsv"
    Map[String, String] ref_index2gsurl = read_map(ref_index_file)
    String genome_url = ref_index2gsurl[genome]

    if (generate_count_config.sample_ids[0] != '') {
        scatter (sample_id in generate_count_config.sample_ids) {
            call mfs.run_merge_fastqs as MergeFastqs {
                input:
                    sample_id = sample_id,
                    fastq_directories = generate_count_config.inpdirs[sample_id],
                    output_directory = output_directory,
                    docker_registry = docker_registry,
                    disk_space = disk_space,
                    zones = zones,
                    memory = memory,
                    preemptible = preemptible
            }

            if (run_count && count_tool == 'StarSolo') {
                call sts.starsolo as StarSolo {
                    input:
                        sample_id = sample_id,
                        r1_fastq = MergeFastqs.fastqs['R1'],
                        r2_fastq = MergeFastqs.fastqs['R2'],
                        genome_url = genome_url + '/starsolo.tar.gz',
                        chemistry = chemistry,
                        output_directory = output_directory,
                        num_cpu = num_cpu,
                        star_version = starsolo_star_version,
                        docker_registry = docker_registry,
                        disk_space = disk_space,
                        preemptible = preemptible,
                        zones = zones,
                        memory = memory
                }
            }

            if (run_count && count_tool == 'Alevin') {
                call ale.alevin as Alevin {
                    input:
                        sample_id = sample_id,
                        r1_fastq = MergeFastqs.fastqs['R1'],
                        r2_fastq = MergeFastqs.fastqs['R2'],
                        genome_url = genome_url + '/alevin.tar.gz',
                        chemistry = chemistry,
                        output_directory = output_directory,
                        num_cpu = num_cpu,
                        docker_registry = docker_registry,
                        alevin_version = alevin_version,
                        disk_space = disk_space,
                        preemptible = preemptible,
                        zones = zones,
                        memory = memory
                }
            }

            if (run_count && count_tool == 'Bustools') {
                call kbc.bustools as Bustools {
                    input:
                        sample_id = sample_id,
                        r1_fastq = MergeFastqs.fastqs['R1'],
                        r2_fastq = MergeFastqs.fastqs['R2'],
                        genome_url = genome_url + '/bustools.tar.gz',
                        chemistry = chemistry,
                        output_directory = output_directory,
                        num_cpu = num_cpu,
                        output_loom = bustools_output_loom,
                        output_h5ad = bustools_output_h5ad,
                        bustools_version = bustools_version,
                        disk_space = disk_space,
                        preemptible = preemptible,
                        zones = zones,
                        memory = memory
                }
            }

            if (run_count && count_tool == 'Optimus') {
                call opm.optimus_count as Optimus {
                    input:
                        sample_id = sample_id,
                        r1_fastq = MergeFastqs.fastqs['R1'],
                        r2_fastq = MergeFastqs.fastqs['R2'],
                        i1_fastq = MergeFastqs.fastqs['I1'],
                        genome_url = genome_url + '/optimus.tar.gz',
                        chemistry = chemistry,
                        output_directory = output_directory,
                        output_loom = optimus_output_loom,
                        docker_registry = docker_registry,
                        version = optimus_version,
                        disk_space = disk_space,
                        preemptible = preemptible,
                        zones = zones,
                        memory = memory
                }
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
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import re, sys
        import pandas as pd

        df = pd.read_csv('${input_tsv_file}', sep = '\t', header = 0, dtype = str, index_col = False)
        for c in df.columns:
            df[c] = df[c].str.strip()

        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['Sample'].str.contains(regex_pat)):
            print('Sample must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)

        with open('sample_ids.txt', 'w') as fo:
            for idx, row in df.iterrows():
                fo.write(row['Sample'] + '\n')

        CODE
    }

    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Map[String, String] inpdirs = read_map(input_tsv_file)
    }

    runtime {
        docker: "${docker_registry}/count"
        zones: zones
        preemptible: "${preemptible}"
    }
}