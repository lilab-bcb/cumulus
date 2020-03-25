version 1.0

workflow souporcell {
    input {
        File input_sample_sheet
        String output_directory
        Int min_num_genes = 500
        String genome
        String souporcell_version = "2020.03"

        String docker_registry = "cumulusprod"
        Int num_cpu = 32
        Int disk_space = 500
        Int memory = 120
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    }

    File ref_index_file = "gs://regev-lab/resources/cellranger/index.tsv"
    # File ref_index_file = "index.tsv"
    Map[String, String] ref_index2gsurl = read_map(ref_index_file)
    String genome_url = ref_index2gsurl[genome]

    call generate_demux_config as DemuxConfig {
        input:
            input_sample_sheet = input_sample_sheet,
            souporcell_version = souporcell_version,
            docker_registry = docker_registry,
            zones = zones,
            preemptible = preemptible
    }

    if (DemuxConfig.sample_ids[0] != '') {
        scatter (sample_id in DemuxConfig.sample_ids) {
            call run_souporcell_demux as SouporcellDemux {
                input:
                    sample_id = sample_id,
                    input_h5 = DemuxConfig.input_h5s[sample_id],
                    input_bam = DemuxConfig.input_bams[sample_id],
                    genome = genome_url,
                    input_vcf_gz = DemuxConfig.input_vcfs[sample_id],
                    min_num_genes = min_num_genes,
                    output_directory = output_directory,
                    donor_names = DemuxConfig.input_donor_names[sample_id],
                    num_clusters = DemuxConfig.input_n_clusters[sample_id],
                    souporcell_version = souporcell_version,
                    docker_registry = docker_registry,
                    num_cpu = num_cpu,
                    disk_space = disk_space,
                    memory = memory,
                    preemptible = preemptible,
                    zones = zones
            }
        }
    }

    output {
        String output_folder = "~{output_directory}"
        Array[File]? demux_output = SouporcellDemux.output_zarr
    }

}

task generate_demux_config {
    input {
        File input_sample_sheet
        String souporcell_version
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

        df = pd.read_csv('~{input_sample_sheet}', header = 0, dtype = str, index_col = False)
        for c in df.columns:
            df[c] = df[c].str.strip()
        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['Sample'].str.contains(regex_pat)):
            print('Sample must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)
            
        with open('run_ids.txt', 'w') as fo_ids, open('h5s.txt', 'w') as fo_h5, open('bams.txt', 'w') as fo_bam, open('vcfs.txt', 'w') as fo_vcf, open('donor_names.txt', 'w') as fo_donors, open('n_clusters.txt', 'w') as fo_clusters:
            for idx, row in df.iterrows():
                fo_ids.write(row['Sample'] + '\n')
                fo_h5.write(row['Sample'] + '\t' + row['H5_file'] + '\n')
                fo_bam.write(row['Sample'] + '\t' + row['Bam_file'] + '\n')
                fo_vcf.write(row['Sample'] + '\t' + row['Genotypes'] + '\n')

                if 'Donors' in df.columns:
                    fo_donors.write(row['Sample'] + '\t' + row['Donors'] + '\n')
                else:
                    fo_donors.write(row['Sample'] + '\tnull\n')
                fo_clusters.write(row['Sample'] + '\t' + row['Clusters'] + '\n')
        CODE
    }

    output {
        Array[String] sample_ids = read_lines('run_ids.txt')
        Map[String, String] input_h5s = read_map('h5s.txt')
        Map[String, String] input_bams = read_map('bams.txt')
        Map[String, String] input_vcfs = read_map('vcfs.txt')
        Map[String, String] input_donor_names = read_map('donor_names.txt')
        Map[String, String] input_n_clusters = read_map('n_clusters.txt')
    }

    runtime {
        docker: "~{docker_registry}/souporcell:~{souporcell_version}"
        zones: zones
        preemptible: "~{preemptible}"
    }
}

task run_souporcell_demux {
    input {
        String sample_id
        File input_h5
        File input_bam
        File genome
        File input_vcf_gz
        Int min_num_genes
        String output_directory
        String donor_names
        Int num_clusters
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
        python /opt/extract_barcodes_for_souporcell.py ~{input_h5} result/~{sample_id}.barcodes.tsv ~{min_num_genes}
        souporcell_pipeline.py -i ~{input_bam} -b result/~{sample_id}.barcodes.tsv -f genome_ref/fasta/genome.fa -t ~{num_cpu} -o result -k ~{num_clusters}

        python <<CODE
        from subprocess import check_call

        call_args = ['python', '/opt/match_donors.py']
        if donor_names is not 'null':
            call_args.extend(['--donor-names', '${donor_names}'])

        call_args.extend(['result/cluster_genotypes.vcf', '~{input_vcf_gz}', 'result/clusters.tsv', '~{input_h5}', 'result/~{sample_id}_demux.zarr'])

        print(' '.join(call_args))
        check_call(call_args, stdout = 'result/match_donors.log')
        CODE

        gsutil -q -m rsync -r result ~{output_directory}/~{sample_id}
        # mkdir -p ~{output_directory}/~{sample_id}
        # cp -r result/* ~{output_directory}/~{sample_id}
    }

    output {
        String output_folder = "~{output_directory}/~{sample_id}"
        File output_zarr = "~{output_directory}/~{sample_id}/~{sample_id}_demux.zarr"
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