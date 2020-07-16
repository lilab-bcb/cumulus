version 1.0

workflow souporcell {
    input {
        String sample_id
        String output_directory
        String input_rna
        String input_bam
        String genome_url
        String ref_genotypes_url
        Boolean de_novo_mode
        Int min_num_genes
        Int num_clusters
        String donor_rename = ''
        String souporcell_version = "2020.06"
        String docker_registry = "cumulusprod"
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
            de_novo_mode = de_novo_mode,
            min_num_genes = min_num_genes,
            num_clusters = num_clusters,
            donor_rename = donor_rename,
            souporcell_version = souporcell_version,
            docker_registry = docker_registry,
            num_cpu = num_cpu,
            disk_space = disk_space,
            memory = memory,
            preemptible = preemptible,
            zones = zones
    }

    output {
        String output_folder = "~{output_directory}"
        File output_zarr = run_souporcell.output_zarr
        File monitoringLog = run_souporcell.monitoringLog
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
        Boolean de_novo_mode
        Int min_num_genes
        Int num_clusters
        String donor_rename
        String souporcell_version

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

        if '~{ref_genotypes}' is not '' and '~{de_novo_mode}' is 'false':
            with open('ref_genotypes.vcf', 'w') as fout:
                check_call(['gunzip', '-k', '~{ref_genotypes}', '-c'], stdout = fout)

            souporcell_call_args.extend(['--known_genotypes', 'ref_genotypes.vcf'])

            assert '~{donor_rename}' is not ''

            name_list = '~{donor_rename}'.split(',')
            souporcell_call_args.extend(['--known_genotypes_sample_names'] + name_list)

        print(' '.join(souporcell_call_args))
        check_call(souporcell_call_args)

        match_call_args = ['python', '/opt/match_donors.py']

        if '~{ref_genotypes}' is not '':
            match_call_args.extend(['--ref-genotypes', '~{ref_genotypes}'])

        if '~{donor_rename}' is not '':
            match_call_args.extend(['--donor-names', '~{donor_rename}'])

        match_call_args.extend(['result/cluster_genotypes.vcf', 'result/clusters.tsv', '~{input_rna}', 'result/~{sample_id}_demux.zarr.zip'])

        print(' '.join(match_call_args))

        with open('match_donors.log', 'w') as fout:
            check_call(match_call_args, stdout = fout)
        CODE

        mkdir buffer
        cp match_donors.log buffer
        cp result/~{sample_id}_demux.zarr.zip buffer
        cp result/clusters.tsv buffer
        cp result/cluster_genotypes.vcf buffer
        gsutil -q -m rsync -r buffer ~{output_directory}/~{sample_id}
        # mkdir -p ~{output_directory}/~{sample_id}
        # cp -r buffer/* ~{output_directory}/~{sample_id}
    }

    output {
        String output_folder = "~{output_directory}/~{sample_id}"
        File output_zarr = "buffer/~{sample_id}_demux.zarr.zip"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/souporcell:~{souporcell_version}"
        zones: zones
        memory: "~{memory}G"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_cpu}"
        preemptible: "~{preemptible}"
    }
}
