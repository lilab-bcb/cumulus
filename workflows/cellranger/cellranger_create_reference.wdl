version 1.0

workflow cellranger_create_reference {
    input {
        # Output directory, gs URL
        String output_directory

        # A sample sheet in CSV format allows users to specify more than 1 genomes to build references (e.g. human and mouse). If a sample sheet is provided, input_fasta, input_gtf, and attributes will be ignored.
        File? input_sample_sheet
        # Input gene annotation file in either GTF or GTF.gz format
        String? input_gtf
        # Input genome reference in either FASTA or FASTA.gz format
        String? input_fasta
        # Genome reference name. New reference will be stored in a folder named genome
        String? genome
        # A list of key:value pairs separated by ;. If this option is not None, cellranger mkgtf will be called to filter the user-provided GTF file. See 10x filter with mkgtf for more details
        String? attributes
        # If we want to build pre-mRNA references, in which we use full length transcripts as exons in the annotation file. We follow 10x build Cell Ranger compatible pre-mRNA Reference Package to build pre-mRNA references
        Boolean pre_mrna = false
        # reference version string
        String? ref_version

        # Which docker registry to use
        String docker_registry = "quay.io/cumulus"
        # 10.0.0, 9.0.1, 8.0.1, 7.2.0
        String cellranger_version = "10.0.0"

        # Disk space in GB
        Int disk_space = 100
        # Google cloud zones
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Number of CPUs
        Int num_cpu = 32
        # Memory string
        String memory = "32G"

        # Number of preemptible tries
        Int preemptible = 2
        # Arn string of AWS queue
        String awsQueueArn = ""
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    # Backend: gcp, aws, local
    Boolean use_gcp = sub(output_directory, "^gs://.+$", "gcp") == "gcp"
    Boolean use_aws = sub(output_directory, "^s3://.+$", "aws") == "aws"
    String backend = (if use_gcp then "gcp"  else (if use_aws then "aws" else "local"))

    call generate_create_reference_config {
        input:
            input_sample_sheet = input_sample_sheet,
            input_gtf_file = input_gtf,
            input_fasta = input_fasta,
            genome = genome,
            attributes = attributes,
            docker_registry = docker_registry,
            cellranger_version = cellranger_version,
            zones = zones,
            preemptible = preemptible,
            awsQueueArn = awsQueueArn,
            backend = backend
    }

    scatter (filt_gtf_row in generate_create_reference_config.filt_gtf_input) {
        call run_filter_gtf {
            input:
                input_gtf_file = filt_gtf_row[0],
                attributes = if length(filt_gtf_row)==2 then filt_gtf_row[1] else '',
                pre_mrna = pre_mrna,
                docker_registry = docker_registry,
                cellranger_version = cellranger_version,
                disk_space = disk_space,
                zones = zones,
                memory = memory,
                preemptible = preemptible,
                awsQueueArn = awsQueueArn,
                backend = backend
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
            preemptible = preemptible,
            awsQueueArn = awsQueueArn,
            backend = backend
    }

    output {
        File output_reference = run_cellranger_mkref.output_reference
    }
}



task generate_create_reference_config {
    input {
        File? input_sample_sheet
        String? input_gtf_file
        String? input_fasta
        String? genome
        String? attributes
        String cellranger_version
        String docker_registry
        String zones
        Int preemptible
        String awsQueueArn
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import pandas as pd

        if '~{input_sample_sheet}' != '':
            df = pd.read_csv('~{input_sample_sheet}', header = 0, dtype = str, index_col = False, keep_default_na = False)
            df.columns = df.columns.str.strip()
            for c in df.columns:
                df[c] = df[c].str.strip()
            if 'Attributes' not in df.columns:
                df['Attributes'] = ''
        else:
            df = pd.DataFrame()
            df['Genome'] = ['~{genome}']
            df['Fasta'] = ['~{input_fasta}']
            df['Genes'] = ['~{input_gtf_file}']
            df['Attributes'] = ['~{attributes}']

        with open('genome_names.txt', 'w') as fo1, open('fasta_files.txt', 'w') as fo2, open('filt_gtf_input.tsv', 'w') as fo3:
            new_genome = []
            for index, row in df.iterrows():
                fo1.write(row['Genome'] + '\n')
                fo2.write(row['Fasta'] + '\n')
                fo3.write(row['Genes'] + '\t' + row['Attributes'] + '\n')
                new_genome.append(row['Genome'])

        concated_genome = '_and_'.join(new_genome)
        if concated_genome == '':
            raise ValueError("Genome attribute must be set!")
        print(concated_genome)
        CODE
    }

    output {
        Array[String] genome_names = read_lines("genome_names.txt")
        Array[String] fasta_files = read_lines("fasta_files.txt")
        Array[Array[String]] filt_gtf_input = read_tsv("filt_gtf_input.tsv")
        String concated_genome = read_string(stdout())
    }

    runtime {
        docker: "~{docker_registry}/cellranger:~{cellranger_version}"
        zones: zones
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}

task run_filter_gtf {
    input {
        File input_gtf_file
        String? attributes
        Boolean pre_mrna

        String docker_registry
        String cellranger_version
        Int disk_space
        String zones
        String memory
        Int preemptible
        String awsQueueArn
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        python <<CODE
        import os
        from subprocess import check_call

        # unzip gtf file if needed.
        input_gtf_file = '~{input_gtf_file}'
        root, ext = os.path.splitext(input_gtf_file)

        if ext == '.gz':
            call_args = ['gunzip', '-f', input_gtf_file]
            print(' '.join(call_args))
            check_call(call_args)
            input_gtf_file = root

        root, ext = os.path.splitext(input_gtf_file)
        file_name = os.path.basename(root)

        output_gtf_file = input_gtf_file # in case no filtering

        if '~{attributes}' != '':
            file_name += '.filt'
            output_gtf_file = file_name + '.gtf'
            call_args = ['cellranger', 'mkgtf', input_gtf_file, output_gtf_file]
            attrs = '~{attributes}'.split(';')
            for attr in attrs:
                call_args.append('--attribute=' + attr)
            print(' '.join(call_args))
            check_call(call_args)
            input_gtf_file = output_gtf_file

        if '~{pre_mrna}' == 'true':
            file_name += '.pre_mrna'
            output_gtf_file = file_name + '.gtf'
            call_args = ['awk', 'BEGIN\\x7BFS="\\\\t"; OFS="\\\\t"\\x7D \\x243 == "transcript" \\x7B\\x243="exon"; print\\x7D', input_gtf_file]
            print(' '.join(call_args) + '> ' + output_gtf_file)
            with open(output_gtf_file, 'w') as fo1:
                check_call(call_args, stdout = fo1)

        call_args = ['mv', output_gtf_file, 'run_filter_gtf_out.gtf']
        print(' '.join(call_args))
        check_call(call_args)
        CODE
    }

    output {
        File output_gtf_file = "run_filter_gtf_out.gtf"
    }

    runtime {
        docker: "~{docker_registry}/cellranger:~{cellranger_version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}

task run_cellranger_mkref {
    input {
        Array[String] genomes
        Array[File] fastas
        Array[File] gtfs
        String output_genome
        String output_dir
        String? ref_version

        String docker_registry
        String cellranger_version
        Int disk_space
        String zones
        Int num_cpu
        String memory
        Int preemptible
        String awsQueueArn
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        python <<CODE
        import os
        from subprocess import check_call

        genome_list = '~{sep="," genomes}'.split(',')
        input_fasta_list = '~{sep="," fastas}'.split(',')
        gtf_list = '~{sep="," gtfs}'.split(',')

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

        mem_digit = ""
        for c in "~{memory}":
            if c.isdigit():
                mem_digit += c
            else:
                break

        call_args.extend(['--nthreads=~{num_cpu}', '--memgb=' + mem_digit, '--localcores=~{num_cpu}', '--localmem=' + mem_digit])
        if '~{ref_version}' != '':
            call_args.append('--ref-version=~{ref_version}')

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        tar -czf ~{output_genome}.tar.gz ~{output_genome}
        strato cp ~{output_genome}.tar.gz "~{output_dir}"/
    }

    output {
        File output_reference = "~{output_genome}.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cellranger:~{cellranger_version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}
