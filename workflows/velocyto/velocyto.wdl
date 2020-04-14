workflow velocyto {
  # The input bam file needs to be sorted by position, this can be achieved running samtools sort mybam.bam -o sorted_bam.bam.
  # cellranger generated bamfiles are already sorted this way.
    File bam
    File sample_name
    File barcodes
    # gs://regev-lab/resources/cellranger/gtf/refdata-cellranger-GRCh38-3.0.0.gtf.gz
    # gs://regev-lab/resources/cellranger/gtf/refdata-cellranger-mm10-1.2.0.gtf.gz
    # gs://regev-lab/resources/cellranger/gtf/refdata-cellranger-mm10-3.0.0.gtf.gz
    File gtf
    # gs://regev-lab/resources/masked_gtf/hg38_rmsk_12_2013.gtf.gz
    # gs://regev-lab/resources/masked_gtf/mm10_rmsk_12_2011.gtf.gz
    File? masked_gtf
    String dtype = "uint16"
    String sort_memory = "15G"
    Int sort_cpu = 2
    Int sort_compression = 4
    Int? cpu = 1
    String? memory = "12G"
    Int? extra_disk_space = 20
    Int? preemptible = 2
    String? version = "0.17"
    String docker_registry = "cumulusprod"
    # docker_registry with trailing slashes stripped
    String docker_registry_stripped = sub(docker_registry, "/+$", "")
    Float sort_disk_multiplier = 3.25
    String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"

    call sort_by_cell_barcode {
        input:
            input_bam=bam,
            sample_name=sample_name,
            extra_disk_space=extra_disk_space,
            num_cpu=sort_cpu,
            memory=sort_memory,
            compression=sort_compression,
            preemptible=preemptible,
            version=version,
            docker_registry=docker_registry,
            sort_disk_multiplier=sort_disk_multiplier,
            zones=zones
    }

    call run_velocyto {
        input:
            sample_name=sample_name,
            filtered_barcodes=barcodes,
            gtf=gtf,
            masked_gtf=masked_gtf,
            cell_sorted_bam=sort_by_cell_barcode.bam,
            position_sorted_bam=bam,
            dtype = dtype,
            num_cpu = cpu,
            memory = memory,
            docker_registry = docker_registry_stripped,
            extra_disk_space = extra_disk_space,
            preemptible = preemptible,
            version = version,
            zones=zones
    }

    output {
       File loom = run_velocyto.loom
    }
}

task sort_by_cell_barcode {
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


    command {
        set -e
        monitor_script.sh > monitoring.log &

        python <<CODE
        import subprocess
        import multiprocessing
        input_bam = '${input_bam}'
        mb_available = int(subprocess.check_output('grep MemAvailable /proc/meminfo'.split()).split()[1]) / 1000
        threads_to_use=multiprocessing.cpu_count()
        mb_to_use = int(mb_available / (1+threads_to_use))
        subprocess.check_call(['samtools', 'sort', '-l', '${compression}', '-m', str(mb_to_use) + 'M', '-@', str(threads_to_use), '-t', 'CB', '-O', 'BAM', '-o', 'cellsorted_${sample_name}.bam', '${input_bam}'])
        CODE
    }

    output {
        File bam = "cellsorted_${sample_name}.bam" #  file must be named cellsorted_name.bam in order to avoid sorting in velocyto
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "${docker_registry}/velocyto:${version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + (ceil(sort_disk_multiplier * size(input_bam, "GB")) + 20) + " HDD"
        cpu: num_cpu
        preemptible: preemptible
    }
}

task run_velocyto {
    # outs/filtered_feature_bc_matrix/barcodes.tsv.gz
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

    command {
        set -e

        monitor_script.sh > monitoring.log &

        mv ${cell_sorted_bam} cellsorted_tmp.bam
        mv ${position_sorted_bam} tmp.bam
        python <<CODE
        
        from subprocess import check_call
        import os
        masked_gtf = '${masked_gtf}'
        gtf = '${gtf}'
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

        call_args = ['velocyto', 'run', '-v', '-t', '${dtype}', '-e', '${sample_name}', '-b', '${filtered_barcodes}', '-o', '.']
        if masked_gtf!='':
            call_args.append('-m')
            call_args.append(masked_gtf)

        # velocyto looks for cellsorted_bam
        call_args.append('tmp.bam')
        call_args.append(gtf)
        check_call(call_args)
        
        CODE
    }

    output {
        File loom = "${sample_name}.loom"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "${docker_registry}/velocyto:${version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(position_sorted_bam, "GB")*2 + size(cell_sorted_bam, "GB")*2 + (2*size(gtf, "GB")) + extra_disk_space)+ " HDD"
        cpu: num_cpu
        preemptible: preemptible
    }
}
