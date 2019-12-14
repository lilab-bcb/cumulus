workflow cellranger_create_reference {
    # Output directory, gs URL
    String output_directory
    File? input_sample_sheet
    String? input_gtf_file
    String? input_fasta
    String? genome
    String? attributes
    Boolean pre_mrna = false
    String? ref_version

    String? docker_registry = "cumulusprod/"
    String? cellranger_version = '3.1.0'
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? num_cpu = 1
    Int? memory = 32

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, “/+$“, “”)

    call generate_create_reference_config {
        input:
            input_sample_sheet = input_sample_sheet,
            input_gtf_file = input_gtf_file,
            input_fasta = input_fasta,
            genome = genome,
            attributes = attributes,
            docker_registry = docker_registry,
            cellranger_version = cellranger_version,
            preemptible = preemptible,
            zones = zones
    }

    scatter (filt_gtf_row in generate_create_reference_config.filt_gtf_input) {
        call run_filter_gtf {
            input:
                input_gtf_file = filt_gtf_row[0],
                attributes = filt_gtf_row[1],
                pre_mrna = pre_mrna,
                docker_registry = docker_registry,
                cellranger_version = cellranger_version,
                disk_space = disk_space,
                zones = zones,
                memory = memory,
                preemptible = preemptible
        }
    }


    call run_cellranger_mkref {
        input:
            genomes = generate_create_reference_config.genome_names,
            fastas = generate_create_reference_config.fasta_files,
            gtfs = run_filter_gtf.output_gtf_file,
            output_genome = generate_create_reference_config.concated_genome,
            output_dir = output_directory_stripped,
            ref_version = ref_version,
            docker_registry = docker_registry,
            cellranger_version = cellranger_version,
            disk_space = disk_space,
            memory = memory,
            num_cpu = num_cpu,
            zones = zones,
            preemptible = preemptible
    }
}



task generate_create_reference_config {
    File? input_sample_sheet
    String? input_gtf_file
    String? input_fasta
    String? genome
    String? attributes
    String cellranger_version
    String docker_registry
    Int preemptible
    String zones

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import pandas as pd

        if '${input_sample_sheet}' is not '':			
            df = pd.read_csv('${input_sample_sheet}', header = 0, dtype = str, index_col=False)
            for c in df.columns:
                df[c] = df[c].str.strip()
        else:
            df = pd.DataFrame()
            df['Genome'] = ['${genome}']
            df['Fasta'] = ['${input_fasta}']
            df['Genes'] = [${input_gtf_file}]
            df['Attributes'] = [${attributes}]
            
        with open('genome_names.txt', 'w') as fo1, open('fasta_files.txt', 'w') as fo2, open('filt_gtf_input.tsv', 'w') as fo3:
            new_genome = []
            for index, row in df.iterrows():
                fo1.write(row['Genome'] + '\n')
                fo2.write(row['Fasta'] + '\n')
                fo3.write(row['Genes'] + '\t' + row['Attributes'] + '\n')
                new_genome.append(row['Genome'])
            print('_and_'.join(new_genome))
        CODE
    }

    output {
        Array[String] genome_names = read_lines("genome_names.txt")
        Array[String] fasta_files = read_lines("fasta_files.txt")
        Array[Array[String]] filt_gtf_input = read_tsv("filt_gtf_input.tsv")
        String concated_genome = read_string(stdout())
    }

    runtime {
        docker: "${docker_registry}cellranger:${cellranger_version}"
        zones: zones
        preemptible: "${preemptible}"
    }
}

task run_filter_gtf {
    File input_gtf_file
    String attributes
    Boolean pre_mrna

    String docker_registry
    String cellranger_version
    Int disk_space
    String zones
    Int memory
    Int preemptible

    String pre_mrna_options = 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}'

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        import os
        from subprocess import check_call

        # unzip gtf file if needed.
        input_gtf_file = '${input_gtf_file}'
        root, ext = os.path.splitext(input_gtf_file)

        if ext == '.gz':
            call_args = ['gunzip', '-f', input_gtf_file]
            print(' '.join(call_args))
            check_call(call_args)
            input_gtf_file = root

        root, ext = os.path.splitext(input_gtf_file)
        file_name = os.path.basename(root)

        output_gtf_file = input_gtf_file # in case no filtering		

        if '${attributes}' is not '':
            file_name += '.filt'
            output_gtf_file = file_name + '.gtf'
            call_args = ['cellranger', 'mkgtf', input_gtf_file, output_gtf_file]
            attrs = '${attributes}'.split(';')
            for attr in attrs:
                call_args.append('--attribute=' + attr)
            print(' '.join(call_args))
            check_call(call_args)
            input_gtf_file = output_gtf_file

        if '${pre_mrna}' is 'true':
            file_name += '.pre_mrna'
            output_gtf_file = file_name + '.gtf'
            call_args = ['awk', '${pre_mrna_options}', input_gtf_file]
            print(' '.join(call_args) + '> ' + output_gtf_file)
            with open(output_gtf_file, 'w') as fo1:
                check_call(call_args, stdout = fo1)

        with open('gtf_file.txt', 'w') as fo2:
            fo2.write(output_gtf_file + '\n')
        CODE
	}

    output {
        File output_gtf_file = read_string('gtf_file.txt')
    }

    runtime {
        docker: "${docker_registry}cellranger:${cellranger_version}"
        zones: zones
        memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: 1
        preemptible: "${preemptible}"
    }
}

task run_cellranger_mkref {
    Array[String] genomes
    Array[File] fastas
    Array[File] gtfs
    String output_genome
    String output_dir
    String ref_version

    String docker_registry
    String cellranger_version
    Int disk_space
    Int preemptible
    String zones
    Int num_cpu
    Int memory

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        import os
        from subprocess import check_call

        genome_list = '${sep="," genomes}'.split(',')
        input_fasta_list = '${sep="," fastas}'.split(',')
        gtf_list = '${sep="," gtfs}'.split(',')

        fasta_list = []
        for fa_file in input_fasta_list:
            root, ext = os.path.splitext(fa_file)
            if ext == '.gz':
                call_args = ['gunzip', '-f', fa_file]
                print(' '.join(call_args))
                check_call(call_args)
                fasta_list.append(root)
            else:
                fasta_list.append(fa_file)

        call_args = ['cellranger', 'mkref']

        for genome, fasta, gtf in zip(genome_list, fasta_list, gtf_list):
            call_args.extend(['--genome=' + genome, '--fasta=' + fasta, '--genes=' + gtf])

        call_args.extend(['--nthreads=${num_cpu}', '--memgb=${memory}'])
        if '${ref_version}' is not '':
            call_args.append('--ref-version=${ref_version}')

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        tar -czf ${output_genome}.tar.gz ${output_genome}
        gsutil cp ${output_genome}.tar.gz ${output_dir}
        # mkdir -p ${output_dir}
        # cp ${output_genome}.tar.gz ${output_dir}
    }

    output {
        String output_reference = '${output_genome}.tar.gz'
        File monitoringLog = 'monitoring.log'
    }

    runtime {
        docker: "${docker_registry}cellranger:${cellranger_version}"
        zones: zones
        memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${preemptible}"
    }
}
