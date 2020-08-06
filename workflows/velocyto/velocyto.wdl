version 1.0


workflow velocyto {
    input {
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        File sample_name
        String sort_memory = "15G"
        String docker_registry = "cumulusprod"
        Int extra_disk_space = 20
        Int preemptible = 2
        String version = "0.17"
        File bam
        Int sort_cpu = 2
        String memory = "12G"
        Int cpu = 1
        File barcodes
        String dtype = "uint16"
        Int sort_compression = 4
        File? masked_gtf
        File gtf
        Float sort_disk_multiplier = 3.25
    }
    String docker_registry_stripped = sub(docker_registry, "/+$", "")
    call run_velocyto {
        input:
            filtered_barcodes = barcodes,
            cell_sorted_bam = sort_by_cell_barcode.bam,
            position_sorted_bam = bam,
            gtf = gtf,
            masked_gtf = masked_gtf,
            dtype = dtype,
            version = version,
            sample_name = sample_name,
            zones = zones,
            num_cpu = cpu,
            memory = memory,
            extra_disk_space = extra_disk_space,
            preemptible = preemptible,
            docker_registry = docker_registry_stripped
    }
    call sort_by_cell_barcode {
        input:
            input_bam = bam,
            sample_name = sample_name,
            extra_disk_space = extra_disk_space,
            zones = zones,
            num_cpu = sort_cpu,
            memory = sort_memory,
            compression = sort_compression,
            preemptible = preemptible,
            version = version,
            docker_registry = docker_registry,
            sort_disk_multiplier = sort_disk_multiplier
    }

    output {
        File loom = run_velocyto.loom
    }
}
task sort_by_cell_barcode {
    input {
        File input_bam
        String sample_name
        Int extra_disk_space
        String zones
        Int num_cpu
        String memory
        Int compression
        Int preemptible
        String version
        String docker_registry
        Float sort_disk_multiplier
    }


    output {
        File bam = "cellsorted_${sample_name}.bam"
        File monitoringLog = "monitoring.log"
    }
    command <<<

        set -e
        monitor_script.sh > monitoring.log &

        python <<CODE
        import subprocess
        import multiprocessing
        input_bam = '~{input_bam}'
        mb_available = int(subprocess.check_output('grep MemAvailable /proc/meminfo'.split()).split()[1]) / 1000
        threads_to_use=multiprocessing.cpu_count()
        mb_to_use = int(mb_available / (1+threads_to_use))
        subprocess.check_call(['samtools', 'sort', '-l', '~{compression}', '-m', str(mb_to_use) + 'M', '-@', str(threads_to_use), '-t', 'CB', '-O', 'BAM', '-o', 'cellsorted_~{sample_name}.bam', '~{input_bam}'])
        CODE

    >>>
    runtime {
        preemptible: preemptible
        bootDiskSizeGb: 12
        disks: "local-disk " + (ceil(sort_disk_multiplier * size(input_bam, "GB")) + 20) + " HDD"
        docker: "${docker_registry}/velocyto:${version}"
        cpu: num_cpu
        zones: zones
        memory: memory
    }

}
task run_velocyto {
    input {
        File filtered_barcodes
        File cell_sorted_bam
        File position_sorted_bam
        File gtf
        File? masked_gtf
        String dtype
        String version
        String sample_name
        String zones
        Int num_cpu
        String memory
        Int extra_disk_space
        Int preemptible
        String docker_registry
    }


    output {
        File loom = "${sample_name}.loom"
        File monitoringLog = "monitoring.log"
    }
    command <<<

        set -e

        monitor_script.sh > monitoring.log &

        mv ~{cell_sorted_bam} cellsorted_tmp.bam
        mv ~{position_sorted_bam} tmp.bam
        python <<CODE

        from subprocess import check_call
        import os
        masked_gtf = '~{masked_gtf}'
        gtf = '~{gtf}'
        if gtf.lower().endswith('.tar.gz'): # assume cellranger genome bundle
        check_call(['mkdir', 'genome_dir'])
        check_call(['tar', 'xf', gtf, '-C', 'genome_dir', '--strip-components', '1'])
        gtf = 'genome_dir/genes/genes.gtf'
        elif gtf.lower().endswith('.gtf.gz'):
        check_call(['gunzip', gtf])
        gtf = gtf[0:-3]
        if masked_gtf.lower().endswith('.gz'):
        check_call(['gunzip', masked_gtf])
        masked_gtf = masked_gtf[0:-3]

        call_args = ['velocyto', 'run', '-v', '-t', '~{dtype}', '-e', '~{sample_name}', '-b', '~{filtered_barcodes}', '-o', '.']
        if masked_gtf!='':
        call_args.append('-m')
        call_args.append(masked_gtf)

        # velocyto looks for cellsorted_bam
        call_args.append('tmp.bam')
        call_args.append(gtf)
        check_call(call_args)

        CODE

    >>>
    runtime {
        preemptible: preemptible
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(position_sorted_bam, "GB") * 2 + size(cell_sorted_bam, "GB") * 2 + (2 * size(gtf, "GB")) + extra_disk_space) + " HDD"
        docker: "${docker_registry}/velocyto:${version}"
        cpu: num_cpu
        zones: zones
        memory: memory
    }

}

