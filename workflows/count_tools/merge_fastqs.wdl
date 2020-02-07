workflow merge_fastqs {
    String input_sample_sheet_tsv
    String output_directory

    String? docker_registry = "cumulusprod"
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? memory = 32

    call set_up_merge_config {
        input:
            input_tsv_file = input_sample_sheet_tsv,
            docker_registry = docker_registry,
            zones = zones,
            preemptible = preemptible
    }

    if (set_up_merge_config.sample_ids[0] != '') {
        scatter (sample_id in set_up_merge_config.sample_ids) {
            call run_merge_fastqs {
                input:
                    sample_id = sample_id,
                    fastq_directories = set_up_merge_config.inpdirs[sample_id],
                    output_directory = output_directory,
                    docker_registry = docker_registry,
                    disk_space = disk_space,
                    zones = zones,
                    memory = memory,
                    preemptible = preemptible
            }
        }
    }

    call generate_merged_count_config {
        input:
            sample_ids = set_up_merge_config.sample_ids,
            read_name_array = run_merge_fastqs.read_name_txt,
            output_directory = output_directory,
            docker_registry = docker_registry,
            zones = zones,
            preemptible = preemptible
    }
 
}

task set_up_merge_config {
    File input_tsv_file
    String docker_registry
    String zones
    Int preemptible

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import os, sys, re
        import pandas as pd

        df = pd.read_csv('${input_tsv_file}', sep = '\t', header = 0, dtype = str, index_col = False)
        for c in df.columns:
            df[c] = df[c].str.strip()

        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['Sample'].str.contains(regex_pat)):
            print('Sample must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)

        with open('sample_ids.txt', 'w') as fo1, open('inpdirs.tsv', 'w') as fo2:
            for idx, row in df.iterrows():
                fo1.write(row['Sample'] + '\n')
                fo2.write(row['Sample'] + '\t' + row['Flowcells'] + '\n')

        CODE
    }

    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Map[String, String] inpdirs = read_map('inpdirs.tsv')
    }

    runtime {
        docker: "${docker_registry}/merge-fastqs"
        zones: zones
        preemptible: "${preemptible}"
    }
}

task run_merge_fastqs {
    String sample_id
    String fastq_directories
    String output_directory
    String docker_registry
    Int disk_space
    String zones
    Int memory
    Int preemptible

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &
        mkdir result

        python <<CODE
        import re, os
        import pandas as pd
        import numpy as np
        from subprocess import check_call

        input_dir_list = list(map(lambda x: x.strip(), "${fastq_directories}".split(',')))
        dir_count = 0
        for directory in input_dir_list:
            directory = re.sub('/+$', '', directory)
            call_args = ['gsutil', '-q', '-m', 'cp', '-r', directory, str(dir_count)]
            # call_args = ['cp', '-r', directory, str(dir_count)]
            print(' '.join(call_args))
            check_call(call_args)

            dir_count += 1

        fastq_files = [f for f in os.listdir(str(0)) if re.match('.*.fastq.gz', f)]
        read_names = pd.Series(list(map(lambda f: f.split('.')[-3].split('_')[-2], fastq_files))).unique()
        with open('read_names.txt', 'w') as fo:
            fo.write('\n'.join(read_names) + '\n')

        fastq_dict = dict()
        for i in range(dir_count):
            fname_pattern_list = os.listdir(str(i))[0].split('.')[0].split('_')
            fname_prefix = '_'.join(fname_pattern_list[:-3]) + '_'
            fname_suffix = '_' + fname_pattern_list[-1] + '.fastq.gz'

            lane_ids = np.sort(pd.Series(list(map(lambda f: f.split('.')[-3].split('_')[-3], os.listdir(str(i))))).unique())

            for rname in read_names:
                f_list = [str(i) + '/' + fname_prefix + lane + '_' + rname + fname_suffix for lane in lane_ids]
                if rname in fastq_dict.keys():
                    fastq_dict[rname].extend(f_list)
                else:
                    fastq_dict[rname] = f_list
            
        for rname in read_names:
            call_args = ['cat']
            call_args.extend(fastq_dict[rname])
            output_fastq = 'result/${sample_id}_' + rname + '.fastq.gz'
            print(' '.join(call_args) + ' > ' + output_fastq)
            with open(output_fastq, 'w') as merge_out:
                check_call(call_args, stdout = merge_out)
        CODE

        gsutil -q -m rsync -r result ${output_directory}/${sample_id}
        gsutil -q -m cp read_names.txt ${output_directory}
        # mkdir -p ${output_directory}/${sample_id}
        # cp result/* ${output_directory}/${sample_id}
        # cp read_names.txt ${output_directory}
    }

    output {
        String read_name_txt = "${output_directory}/read_names.txt"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "${docker_registry}/merge-fastqs"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ${disk_space} HDD"
        cpu: 1
        preemptible: "${preemptible}"
    }
}

task generate_merged_count_config {
    Array[String] sample_ids
    Array[String] read_name_array
    String output_directory
    String docker_registry
    String zones
    Int preemptible

    String read_name_txt = read_name_array[0]

    command {
        set -e
        export TMPDIR=/tmp

        gsutil -q -m cp ${read_name_txt} read_names.txt
        # cp ${read_name_txt} read_names.txt

        python <<CODE
        import pandas as pd

        sample_names = '${sep="," sample_ids}'.split(',')

        df = pd.read_csv('read_names.txt', header = None)
        read_names = df[0].values

        with open('count_matrix.csv', 'w') as fo:
            fo.write('Sample,' + ','.join([rname + '_fastq' for rname in read_names]) + '\n')
            for sample in sample_names:
                fo.write(sample + ',' + ','.join(['${output_directory}/sample/sample_' + rname + '.fastq.gz' for rname in read_names]) + '\n')
        CODE

        gsutil -q -m cp count_matrix.csv ${output_directory}
        # cp count_matrix.csv ${output_directory}
    }

    output {
        String count_matrix = "${output_directory}/count_matrix.csv"
    }

    runtime {
        docker: "${docker_registry}/merge-fastqs"
        zones: zones
        preemptible: "${preemptible}"
    }

}