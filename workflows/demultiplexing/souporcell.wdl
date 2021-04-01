version 1.0

workflow souporcell {
    input {
        String sample_id
        String output_directory
        String input_rna
        String input_bam
        String genome_url
        String ref_genotypes_url
        File? common_variants
        Boolean skip_remap
        Boolean de_novo_mode
        Int min_num_genes
        Int num_clusters
        String donor_rename = ''
        String souporcell_version = "2020.07"
        String docker_registry = "quay.io/cumulus"
        Int num_cpu = 32
        Int disk_space = 500
        Int memory = 120
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    }

    if (ref_genotypes_url != 'null') {
        File ref_genotypes = ref_genotypes_url
    }

    call run_souporcell {
        input:
            sample_id = sample_id,
            output_directory = output_directory,
            input_rna = input_rna,
            input_bam = input_bam,
            genome = genome_url,
            ref_genotypes = ref_genotypes,
            common_variants = common_variants,
            skip_remap = skip_remap,
            donor_rename = donor_rename,
            de_novo_mode = de_novo_mode,
            min_num_genes = min_num_genes,
            num_clusters = num_clusters,
            version = souporcell_version,
            docker_registry = docker_registry,
            num_cpu = num_cpu,
            disk_space = disk_space,
            memory = memory,
            preemptible = preemptible,
            zones = zones
    }

    call match_donors {
        input:
            sample_id = sample_id,
            output_directory = output_directory,
            input_rna = input_rna,
            souporcell_cluster_tsv = run_souporcell.souporcell_cluster_tsv,
            souporcell_genotypes_vcf = run_souporcell.souporcell_genotypes_vcf,
            ref_genotypes = ref_genotypes,
            donor_rename = donor_rename,
            docker_registry = docker_registry,
            version = souporcell_version,
            zones = zones,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible
    }

    output {
        String output_folder = match_donors.output_folder
        File output_zarr = match_donors.output_zarr
        Array[File] monitoringLog = [run_souporcell.monitoringLog, match_donors.monitoringLog]
    }

}

task run_souporcell {
    input {
        String sample_id
        String output_directory
        File input_rna
        File input_bam
        File genome
        File? ref_genotypes
        File? common_variants
        Boolean skip_remap
        String? donor_rename
        Boolean de_novo_mode
        Int min_num_genes
        Int num_clusters
        String version

        String docker_registry
        Int num_cpu
        Int disk_space
        Int memory
        Int preemptible
        String zones
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        mkdir genome_ref
        tar -zxf "~{genome}" -C genome_ref --strip-components 1
        rm "~{genome}"

        mkdir result
        python /opt/extract_barcodes_from_rna.py ~{input_rna} result/~{sample_id}.barcodes.tsv ~{min_num_genes}

        python <<CODE
        from subprocess import check_call

        souporcell_call_args = ['souporcell_pipeline.py', '-i', '~{input_bam}', '-b', 'result/~{sample_id}.barcodes.tsv', '-f', 'genome_ref/fasta/genome.fa', '-t', '~{num_cpu}', '-o', 'result', '-k', '~{num_clusters}']

        if '~{common_variants}' is not '':
            file_ext = '~{common_variants}'.split('.')[-1]
            if file_ext == 'gz':
                with open('common_variants.vcf', 'w') as fout:
                    check_call(['gunzip', '~{common_variants}', '-c'], stdout = fout)
            else:
                check_call(['mv', '~{common_variants}', 'common_variants.vcf'])
            souporcell_call_args.extend(['--common_variants', 'common_variants.vcf'])

            if '~{skip_remap}' is 'true':
                souporcell_call_args.append('--skip_remap')

        if '~{ref_genotypes}' is not '' and '~{de_novo_mode}' is 'false':
            file_ext = '~{ref_genotypes}'.split('.')[-1]
            if file_ext == 'gz':
                with open('ref_genotypes.vcf', 'w') as fout:
                    check_call(['gunzip', '-k', '~{ref_genotypes}', '-c'], stdout = fout)
            else:
                check_call(['mv', '~{ref_genotypes}', 'ref_genotypes.vcf'])

            souporcell_call_args.extend(['--known_genotypes', 'ref_genotypes.vcf'])

            assert '~{donor_rename}' is not ''

            name_list = '~{donor_rename}'.split(',')
            souporcell_call_args.extend(['--known_genotypes_sample_names'] + name_list)

        print(' '.join(souporcell_call_args))
        check_call(souporcell_call_args)

        CODE
    }

    output {
        File souporcell_cluster_tsv = "result/clusters.tsv"
        File souporcell_genotypes_vcf = "result/cluster_genotypes.vcf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/souporcell:~{version}"
        zones: zones
        memory: "~{memory}G"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_cpu}"
        preemptible: "~{preemptible}"
    }
}

task match_donors {
    input {
        String sample_id
        String output_directory
        File input_rna
        File souporcell_cluster_tsv
        File souporcell_genotypes_vcf
        File? ref_genotypes
        String? donor_rename

        String docker_registry
        String version
        String zones
        Int memory
        Int disk_space
        Int preemptible
    }

    Float file_size = size([input_rna, souporcell_cluster_tsv, souporcell_genotypes_vcf], "GB")

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        from subprocess import check_call

        call_args = ['python', '/opt/match_donors.py']

        if '~{ref_genotypes}' is not '':
            call_args.extend(['--ref-genotypes', '~{ref_genotypes}'])

        if '~{donor_rename}' is not '':
            call_args.extend(['--donor-names', '~{donor_rename}'])

        call_args.extend(['~{souporcell_genotypes_vcf}', '~{souporcell_cluster_tsv}', '~{input_rna}', '~{sample_id}_demux.zarr.zip'])

        print(' '.join(call_args))

        with open('match_donors.log', 'w') as fout:
            check_call(call_args, stdout = fout)

        CODE

        mkdir result
        mv match_donors.log ~{sample_id}_demux.zarr.zip ~{souporcell_cluster_tsv} ~{souporcell_genotypes_vcf} result/
        gsutil -q -m rsync -r result ~{output_directory}/~{sample_id}
        # mkdir -p ~{output_directory}/~{sample_id}
        # cp -r result/* ~{output_directory}/~{sample_id}
    }

    output {
        String output_folder = "~{output_directory}/~{sample_id}"
        File output_zarr = "result/~{sample_id}_demux.zarr.zip"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/souporcell:~{version}"
        zones: zones
        memory: "~{memory}G"
        disks: "local-disk ~{disk_space} HDD"
        preemptible: "~{preemptible}"
    }
}
