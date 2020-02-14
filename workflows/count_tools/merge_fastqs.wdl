version 1.0

task run_merge_fastqs {
    input {
        String sample_id
        String fastq_directories
        String output_directory
        String docker_registry
        Int disk_space
        String zones
        Int memory
        Int preemptible
    } 

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
            call_args = ['mkdir', '-p', str(dir_count)]
            print(' '.join(call_args))
            check_call(call_args)

            directory = re.sub('/+$', '', directory)
            call_args = ['gsutil', '-q', '-m', 'rsync', '-r', directory, str(dir_count)]
            # call_args = ['rsync', '-r', directory + '/', str(dir_count)]
            print(' '.join(call_args))
            check_call(call_args)

            dir_count += 1

        fastq_files = [f for f in os.listdir(str(0)) if re.match('.*.fastq.gz', f)]
        read_names = pd.Series(list(map(lambda f: f.split('.')[-3].split('_')[-2], fastq_files))).unique()
        with open('fastqs.tsv', 'w') as fo:
            for rname in read_names:
                fastq_file_name = '${output_directory}/merged_fastqs/${sample_id}/${sample_id}_' + rname + '.fastq.gz'
                fo.write(rname + '\t' + fastq_file_name + '\n')

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

        gsutil -q -m rsync -r result ${output_directory}/merged_fastqs/${sample_id}
        # mkdir -p ${output_directory}/merged_fastqs/${sample_id}
        # cp result/* ${output_directory}/merged_fastqs/${sample_id}
    }

    output {
        Map[String, String] fastqs = read_map("fastqs.tsv")
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "${docker_registry}/count"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ${disk_space} HDD"
        cpu: 3
        preemptible: "${preemptible}"
    }
}