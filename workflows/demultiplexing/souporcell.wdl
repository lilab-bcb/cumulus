version 1.0

workflow souporcell {
    input {
        # Sample ID
        String sample_id
        # Output directory (gs url + path)
        String output_directory
        # Link to a RNA count matrix
        String input_rna
        # Link to a position-sorted BAM
        String input_bam
        # Link to a reference genome
        String genome_url
        # Link to reference genotype files
        String ref_genotypes_url
        # Number of expected clusters when doing clustering
        Int num_clusters
        # If true, run souporcell in de novo mode without reference genotypes
        Boolean de_novo_mode
        # Users can provide a common variants list in VCF format for Souporcell to use, instead of calling SNPs de novo
        File? common_variants
        # Skip remap step. Only recommended in non denovo mode or common variants are provided
        Boolean skip_remap
        # Set if your umi tag is not UB
        String umi_tag
        # A comma-separated list of donor names for renaming clusters achieved by souporcell
        String donor_rename
        # Only demultiplex cells/nuclei with at least <min_num_genes> expressed genes
        Int min_num_genes
        # Souporcell version to use
        String souporcell_version
        # Which docker registry to use
        String docker_registry
        # Google cloud zones
        String zones
        # Number of CPUs to request for souporcell per pair
        Int num_cpu
        # Disk space (integer) in GB needed for souporcell per pair
        Int disk_space
        # Memory size string for souporcell per pair
        String memory
        # Number of preemptible tries
        Int preemptible
        # Arn string of AWS queue
        String awsQueueArn
        # Backend
        String backend
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
            umi_tag = umi_tag,
            de_novo_mode = de_novo_mode,
            min_num_genes = min_num_genes,
            num_clusters = num_clusters,
            version = souporcell_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            disk_space = disk_space,
            memory = memory,
            preemptible = preemptible,
            awsQueueArn = awsQueueArn,
            backend = backend
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
            preemptible = preemptible,
            awsQueueArn = awsQueueArn,
            backend = backend
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
        String umi_tag
        Boolean de_novo_mode
        Int min_num_genes
        Int num_clusters

        String version
        String docker_registry
        String zones
        Int num_cpu
        Int disk_space
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

        mkdir genome_ref
        tar -zxf "~{genome}" -C genome_ref --strip-components 1
        rm "~{genome}"

        mkdir result
        python /opt/extract_barcodes_from_rna.py ~{input_rna} result/~{sample_id}.barcodes.tsv ~{min_num_genes}

        python <<CODE
        from subprocess import check_call

        souporcell_call_args = ['souporcell_pipeline.py', '-i', '~{input_bam}', '-b', 'result/~{sample_id}.barcodes.tsv', '-f', 'genome_ref/fasta/genome.fa', '-t', '~{num_cpu}', '-o', 'result', '-k', '~{num_clusters}']

        if '~{de_novo_mode}' == 'false':
            assert '~{ref_genotypes}' != '', "Reference mode must have a reference genotype vcf file provided!"
            file_ext = '~{ref_genotypes}'.split('.')[-1]
            if file_ext == 'gz':
                with open('ref_genotypes.vcf', 'w') as fout:
                    check_call(['gunzip', '-k', '~{ref_genotypes}', '-c'], stdout = fout)
            else:
                check_call(['mv', '~{ref_genotypes}', 'ref_genotypes.vcf'])

            souporcell_call_args.extend(['--known_genotypes', 'ref_genotypes.vcf'])
        else:
            assert '~{de_novo_mode}' == 'true'
            if '~{common_variants}' != '':
                file_ext = '~{common_variants}'.split('.')[-1]
                if file_ext == 'gz':
                    with open('common_variants.vcf', 'w') as fout:
                        check_call(['gunzip', '~{common_variants}', '-c'], stdout = fout)
                else:
                    check_call(['mv', '~{common_variants}', 'common_variants.vcf'])
                souporcell_call_args.extend(['--common_variants', 'common_variants.vcf'])

        if '~{skip_remap}' == 'true':
            if '~{de_novo_mode}' == 'true' and '~{common_variants}' == '':
                print("Warning: if de novo mode is true and no common variants provided, skip remap is not recommended and thus is turned off!")
            else:
                check_call(['samtools', 'index', '-@', '~{num_cpu}', '~{input_bam}'])
                souporcell_call_args.extend(['--skip_remap', 'True'])

        if '~{umi_tag}' not in ["UB", ""]:
            souporcell_call_args.extend(["--umi_tag", "~{umi_tag}"])

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
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        queueArn: awsQueueArn
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

        String version
        String docker_registry
        String zones
        String memory
        Int disk_space
        Int preemptible
        String awsQueueArn
        String backend
    }

    Float file_size = size([input_rna, souporcell_cluster_tsv, souporcell_genotypes_vcf], "GB")

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        python <<CODE
        from subprocess import check_call

        call_args = ['python', '/opt/match_donors.py']

        if '~{ref_genotypes}' != '':
            call_args.extend(['--ref-genotypes', '~{ref_genotypes}'])

        if '~{donor_rename}' != '':
            call_args.extend(['--donor-names', '~{donor_rename}'])

        call_args.extend(['~{souporcell_genotypes_vcf}', '~{souporcell_cluster_tsv}', '~{input_rna}', '~{sample_id}_demux.zarr.zip'])

        print(' '.join(call_args))

        with open('match_donors.log', 'w') as fout:
            check_call(call_args, stdout = fout)

        CODE

        mkdir result
        mv match_donors.log "~{sample_id}"_demux.zarr.zip ~{souporcell_cluster_tsv} ~{souporcell_genotypes_vcf} result/
        strato sync -m result "~{output_directory}/~{sample_id}"
    }

    output {
        String output_folder = "~{output_directory}/~{sample_id}"
        File output_zarr = "result/~{sample_id}_demux.zarr.zip"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/souporcell:~{version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}
