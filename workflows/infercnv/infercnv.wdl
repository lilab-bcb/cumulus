version 1.0

workflow InferCNV {
    input {
        File input_zarr
        String sample_id
        String output_directory
        File gene_ordering_csv
        String protocol = "tenX"
        Boolean cluster_by_groups = true
        Boolean HMM = false
        String? ref_group_names

        String docker_registry = "cumulusbeta"
        String inferCNV_version = "1.5.1"
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        Int memory = 20
        Int disk_space = 50
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

    call preprocess {
        input:
            input_zarr = input_zarr,
            sample_id = sample_id,
            docker_registry = docker_registry,
            version = inferCNV_version,
            zones = zones,
            preemptible = preemptible,
            memory = memory,
            disk_space = disk_space
    }

    call run_inferCNV {
        input:
            sample_id = sample_id,
            output_directory = output_directory_stripped,
            input_raw_count = preprocess.raw_count,
            gene_ordering_csv = gene_ordering_csv,
            sample_annotation_csv = preprocess.sample_annotation,
            protocol = protocol,
            cluster_by_groups = cluster_by_groups,
            HMM = HMM,
            ref_group_names = ref_group_names,
            docker_registry = docker_registry,
            version = inferCNV_version,
            zones = zones,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible
    }

    call gen_CNV_metrics {
        input:
            sample_id = sample_id,
            output_directory = output_directory_stripped,
            input_cnv_obj = run_inferCNV.cnv_obj,
            anno_zarr = input_zarr,
            sample_annotation_csv = preprocess.sample_annotation,
            docker_registry = docker_registry,
            version = inferCNV_version,
            zones = zones,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible
    }


    output {
        String output_folder = run_inferCNV.output_cnv_folder
        Array[File] output_cnv_plots = flatten([run_inferCNV.output_cnv_plots, gen_CNV_metrics.output_plots])
        File output_raw_count = preprocess.raw_count
        File output_cnv_zarr = gen_CNV_metrics.output_cnv_zarr
    }
}

task preprocess {
    input {
        File input_zarr
        String sample_id
        String docker_registry
        String version
        String zones
        Int preemptible
        Int memory
        Int disk_space
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python /software/preprocess.py ~{input_zarr} ~{sample_id}
        gzip ~{sample_id}.csv

    }

    output {
        File raw_count = '~{sample_id}.csv.gz'
        File sample_annotation = '~{sample_id}.sample_annotation.csv'
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/infercnv:~{version}"
        zones: zones
        preemptible: preemptible
        memory: "~{memory}G"
        disks: "local-disk ~{disk_space} HDD"
    }
}

task run_inferCNV {
    input {
        String sample_id
        String output_directory
        File input_raw_count
        File gene_ordering_csv
        File sample_annotation_csv
        String protocol
        Boolean cluster_by_groups
        Boolean HMM
        String? ref_group_names

        String docker_registry
        String version
        String zones
        Int memory
        Int disk_space
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        from subprocess import check_call

        call_args = ['Rscript', '/software/run_inferCNV.R', '--out-dir=inferCNV_result', '--raw-counts-matrix=~{input_raw_count}', '--gene-order-file=~{gene_ordering_csv}', '--annotations-file=~{sample_annotation_csv}']
        if '~{protocol}' is not '':
            call_args.append('--protocol=~{protocol}')
        if '~{cluster_by_groups}' is 'true':
            call_args.append('--cluster-by-groups')
        if '~{HMM}' is 'true':
            call_args.append('--HMM')
        if '~{ref_group_names}' is not '':
            call_args.append("--ref-group-names=~{ref_group_names}")

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        gsutil -q -m rsync -r inferCNV_result ~{output_directory}/~{sample_id}/inferCNV_result
        # mkdir -p ~{output_directory}/~{sample_id}
        # cp -r inferCNV_result ~{output_directory}/~{sample_id}/
    }

    output {
        String output_cnv_folder = '~{output_directory}/~{sample_id}/inferCNV_result'
        File cnv_obj = 'inferCNV_result/12_remove_ref_avg_from_obs_adjust.infercnv_obj'
        Array[File] output_cnv_plots = ['inferCNV_result/infercnv.preliminary.png', 'inferCNV_result/infercnv.png']
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/infercnv:~{version}"
        zones: zones
        memory: "~{memory}G"
        disks: "local-disk ~{disk_space} HDD"
        preemptible: preemptible
    }
}

task gen_CNV_metrics {
    input {
        String sample_id
        String output_directory
        File input_cnv_obj
        File anno_zarr
        File sample_annotation_csv
        String docker_registry
        String version
        String zones
        Int memory
        Int disk_space
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        Rscript /software/calc_cnv_metrics.R ~{input_cnv_obj} ~{sample_id}_cnv
        python /software/insert_cnv_scores.py ~{anno_zarr} ~{sample_id}_cnv.csv ~{sample_annotation_csv} ~{sample_id}
        pegasus plot scatter --basis umap --attributes infercnv_cell_types --wspace 1.6 ~{sample_id}.cnv.zarr.zip ~{sample_id}.ct.umap.pdf
        pegasus plot scatter --basis umap --attributes sqrt_cnv ~{sample_id}.cnv.zarr.zip ~{sample_id}.cnv_scores.umap.pdf

        python <<CODE
        import pegasusio as io
        import matplotlib.pyplot as plt
        import seaborn as sns

        data = io.read_input("~{sample_id}.cnv.zarr.zip")
        plt.figure(figsize = (8, 6))
        plt.subplots_adjust(hspace = 0.3, bottom = 0.3)
        plt.title("Square-root of CNV products of sample ~{sample_id}")
        ax = sns.violinplot(x = 'infercnv_cell_types', y = 'sqrt_cnv', data = data.obs, cut = 0, linewidth = 1)
        ax.set_ylabel("Sqrt CNV product")
        ax.set_xlabel("Sample-specific cell types")
        ax.set_xticklabels(ax.get_xticklabels(), rotation = -90)
        plt.savefig("~{sample_id}.cnv.pdf", dpi = 300)
        plt.close()
        CODE

        mkdir result
        cp ~{sample_id}.cnv.zarr.zip *.pdf result/
        gsutil -q -m rsync -r result ~{output_directory}/~{sample_id}/cnv_result
        # mkdir -p ~{output_directory}/~{sample_id}/cnv_result
        # cp result/* ~{output_directory}/~{sample_id}/cnv_result/
    }

    output {
        File output_cnv_zarr = "~{sample_id}.cnv.zarr.zip"
        Array[File] output_plots = ["~{sample_id}.ct.umap.pdf", "~{sample_id}.cnv_scores.umap.pdf", "~{sample_id}.cnv.pdf"]
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/infercnv:~{version}"
        zones: zones
        memory: "~{memory}G"
        disks: "local-disk ~{disk_space} HDD"
        preemptible: preemptible
    }
}
