version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:count_tools_merge_fastqs/versions/5/plain-WDL/descriptor" as mfs
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:count_tools_starsolo/versions/10/plain-WDL/descriptor" as sts
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:count_tools_bustools/versions/4/plain-WDL/descriptor" as kbc
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:count_tools_alevin/versions/2/plain-WDL/descriptor" as ale
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:count_tools_optimus/versions/7/plain-WDL/descriptor" as opm

#import "merge_fastqs.wdl" as mfs
#import "star-solo.wdl" as sts
#import "bustools.wdl" as kbc
#import "alevin.wdl" as ale
#import "optimus-count.wdl" as opm

workflow count {
    input {
        File input_tsv_file

        String genome
        String chemistry
        String output_directory

        # Count
        Boolean run_count = true
        String count_tool = "StarSolo"

        String docker_registry = "cumulusprod"
        Int num_cpu = 32
        Int disk_space = 500
        Int memory = 120
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"

        # Merge-Fastqs
        Int merge_fastqs_memory = 32

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
        Boolean optimus_output_loom = true
    }

    call generate_count_config as Config {
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

    if (Config.sample_ids[0] != '') {
        scatter (sample_id in Config.sample_ids) {
            Boolean do_merge = if Config.merge_flags[sample_id] == 'true' then true else false

            if (do_merge) {
                call mfs.run_merge_fastqs as MergeFastqs {
                    input:
                        sample_id = sample_id,
                        fastq_directories = Config.inpdirs[sample_id],
                        output_directory = output_directory,
                        docker_registry = docker_registry,
                        disk_space = disk_space,
                        zones = zones,
                        memory = merge_fastqs_memory,
                        preemptible = preemptible
                }
            }

            Map[String, String] fastqs = select_first([MergeFastqs.fastqs, Config.samples_no_merge[sample_id]])

            if (run_count && count_tool == 'StarSolo') {
                call sts.starsolo as StarSolo {
                    input:
                        sample_id = sample_id,
                        r1_fastq = fastqs['R1'],
                        r2_fastq = fastqs['R2'],
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
                        r1_fastq = fastqs['R1'],
                        r2_fastq = fastqs['R2'],
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
                        r1_fastq = fastqs['R1'],
                        r2_fastq = fastqs['R2'],
                        genome_url = genome_url + '/bustools.tar.gz',
                        chemistry = chemistry,
                        output_directory = output_directory,
                        num_cpu = num_cpu,
                        output_loom = bustools_output_loom,
                        output_h5ad = bustools_output_h5ad,
                        bustools_version = bustools_version,
                        docker_registry = docker_registry,
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
                        r1_fastq = fastqs['R1'],
                        r2_fastq = fastqs['R2'],
                        i1_fastq = fastqs['I1'],
                        genome_url = genome_url + '/optimus.tar.gz',
                        chemistry = chemistry,
                        output_directory = output_directory,
                        output_loom = optimus_output_loom,
                        docker_registry = docker_registry,
                        disk_space = disk_space,
                        preemptible = preemptible,
                        cpu_for_copy = num_cpu,
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
        import json, os, re, sys
        import pandas as pd
        from subprocess import check_call

        df = pd.read_csv('~{input_tsv_file}', sep = '\t', header = 0, dtype = str, index_col = False)
        for c in df.columns:
            df[c] = df[c].str.strip()

        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['Sample'].str.contains(regex_pat)):
            print('Sample must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)

        data_no_merge = dict()
        with open('sample_ids.txt', 'w') as fo1, open('inpdirs.tsv', 'w') as fo2, open('merge_flags.tsv', 'w') as fo3:
            data_no_merge = dict()
            for idx, row in df.iterrows():
                fo1.write(row['Sample'] + '\n')
                fo2.write(row['Sample'] + '\t' + row['Flowcells'] + '\n')

                input_dir_list = list(map(lambda x: x.strip(), row['Flowcells'].split(',')))
                if len(input_dir_list) > 1:
                    fo3.write(row['Sample'] + '\t' + "true\n")
                    read_item = dict()
                    read_item['read'] = "merged"
                    data_no_merge[row['Sample']] = read_item
                else:
                    dir_name = input_dir_list[0] if input_dir_list[0][-1] != '/' else input_dir_list[0][:-1]
                    with open(row['Sample'] + "_tmp.txt", 'w') as tmp_fout:
                        check_call(['gsutil', 'ls', dir_name], stdout = tmp_fout)
                        # check_call(['ls', dir_name], stdout = tmp_fout)
                    with open(row['Sample'] + "_tmp.txt", 'r') as tmp_fin:
                        file_list = [os.path.basename(line.rstrip('\n')) for line in tmp_fin.readlines()]
                        file_list = [f for f in file_list if f] # Remove empty strings.
                        read_names = pd.Series(list(map(lambda f: f.split('.')[-3].split('_')[-2], file_list))).unique()
                        if read_names.size == len(file_list):
                            # No need to merge fastqs
                            fo3.write(row['Sample'] + '\t' + "false\n")
                            read_item = dict()
                            for rname in read_names:
                                read_item[rname] = dir_name + '/' + [f for f in file_list if re.match('.*_' + rname + '_[^_]+.fastq.gz', f)][0]
                            data_no_merge[row['Sample']] = read_item
                        else:
                            fo3.write(row['Sample'] + '\t' + "true\n")
                            read_item = dict()
                            read_item['read'] = "merged"
                            data_no_merge[row['Sample']] = read_item

        with open('samples_no_merge.json', 'w') as json_fo:
            json.dump(data_no_merge, json_fo)

        CODE
    }

    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Map[String, String] inpdirs = read_map('inpdirs.tsv')
        Map[String, String] merge_flags = read_map('merge_flags.tsv')
        Map[String, Map[String, String]] samples_no_merge = read_json('samples_no_merge.json')
    }

    runtime {
        docker: "~{docker_registry}/count"
        zones: zones
        preemptible: "~{preemptible}"
    }
}
