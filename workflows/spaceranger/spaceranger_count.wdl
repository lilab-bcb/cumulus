version 1.0

workflow spaceranger_count {
    input {
        # Sample ID
        String sample_id
        # A comma-separated list of input FASTQs directories (gs urls)
        String input_fastqs_directories
        # spaceranger output directory, gs url
        String output_directory

        # A reference genome name or a URL to a tar.gz file
        String genome
        # Probe set for FFPE samples, choosing from human_probe_v1, mouse_probe_v1 or a user-provided csv file. Default to '', not FFPE
        String probeset = ""

        # Referece index TSV
        File acronym_file
        # An empty full representing null
        File null_file

        # Brightfield tissue H&E image in .jpg or .tiff format.
        File? image
        # Multi-channel, dark-background fluorescence image as either a single, multi-layer .tiff file, multiple .tiff or .jpg files, or a pre-combined color .tiff or .jpg file.
        Array[File]? darkimage
        # A semi-colon ';' separated string denoting all dark images. This option is equivalent to darkimage and should only be used by spaceranger_workflow
        String? darkimagestr
        #A color composite of one or more fluorescence image channels saved as a single-page, single-file color .tiff or .jpg.
        File? colorizedimage
        # Visium slide serial number.
        String? slide
        # Visium capture area identifier. Options for Visium are A1, B1, C1, D1.
        String? area
        # Slide layout file indicating capture spot and fiducial spot positions.
        File? slidefile
        # Use with automatic image alignment to specify that images may not be in canonical orientation with the hourglass in the top left corner of the image. The automatic fiducial alignment will attempt to align any rotation or mirroring of the image.
        Boolean reorient_images = false
        # Alignment file produced by the manual Loupe alignment step. A --image must be supplied in this case.
        File? loupe_alignment

        # Target panel CSV for targeted gene expression analysis
        File? target_panel

        # If generate bam outputs
        Boolean no_bam = false
        # Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
        Boolean secondary = false

        # spaceranger version
        String spaceranger_version
        # Which docker registry to use: cumulusprod (default) or quay.io/cumulus
        String docker_registry

        # Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per spaceranger job
        Int num_cpu = 32
        # Memory string, e.g. 120G
        String memory = "120G"
        # Disk space in GB
        Int disk_space = 500
        # Number of preemptible tries
        Int preemptible = 2
        # Number of maximum retries when running on AWS
        Int awsMaxRetries = 5
        # Backend
        String backend
    }

    Map[String, String] acronym2url = read_map(acronym_file)

    # If reference is a url
    Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"
    # Replace name with actual url
    File genome_file = (if is_url then genome else acronym2url[genome])

    # If replace probset with its corresponding URL for CSV file
    File probe_file = (if probeset == "" then null_file else (if sub(probeset, "^.+\\.csv$", "CSV") != "CSV" then acronym2url[probeset] else probeset))


    call run_spaceranger_count {
        input:
            sample_id = sample_id,
            input_fastqs_directories = input_fastqs_directories,
            output_directory = output_directory,
            genome_file = genome_file,
            probe_file = probe_file,
            image = image,
            darkimage = darkimage,
            darkimagestr = darkimagestr,
            colorizedimage = colorizedimage,
            slide = slide,
            area = area,
            slidefile = slidefile,
            reorient_images = reorient_images,
            loupe_alignment = loupe_alignment,
            target_panel = target_panel,
            no_bam = no_bam,
            secondary = secondary,
            spaceranger_version = spaceranger_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend
    }

    output {
        String output_count_directory = run_spaceranger_count.output_count_directory
        String output_metrics_summary = run_spaceranger_count.output_metrics_summary
        String output_web_summary = run_spaceranger_count.output_web_summary
        File monitoringLog = run_spaceranger_count.monitoringLog
    }
}

task run_spaceranger_count {
    input {
        String sample_id
        String input_fastqs_directories
        String output_directory
        File genome_file
        File probe_file
        File? image
        Array[File]? darkimage
        String? darkimagestr
        File? colorizedimage
        String? slide
        String? area
        File? slidefile
        Boolean reorient_images
        File? loupe_alignment
        File? target_panel
        Boolean no_bam
        Boolean secondary
        String spaceranger_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &
        mkdir -p genome_dir
        tar xf ~{genome_file} -C genome_dir --strip-components 1

        python <<CODE
        import os
        import re
        import sys
        from subprocess import check_call

        fastqs = []
        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = re.sub('/+$', '', directory) # remove trailing slashes
            call_args = ['strato', 'sync', '--backend', '~{backend}', '-m', directory + '/~{sample_id}', '~{sample_id}_' + str(i)]
            print(' '.join(call_args))
            check_call(call_args)
            fastqs.append('~{sample_id}_' + str(i))

        call_args = ['spaceranger', 'count', '--id=results', '--transcriptome=genome_dir', '--fastqs=' + ','.join(fastqs), '--sample=~{sample_id}', '--jobmode=local']

        def not_null(input_file):
            return (input_file is not '') and (os.path.basename(input_file) != 'null')

        def get_darkimages(darkimage, darkimagestr):
            darkimages = []
            if darkimage != '':
                darkimages = darkimage.split(';')
            elif darkimagestr != '':
                for i, file in enumerate(darkimagestr.split(';')):
                    local_file = '_' + str(i) + '_' + os.path.basename(file)
                    call_args = ['strato', 'cp', '--backend', '~{backend}', file, local_file]
                    print(' '.join(call_args))
                    check_call(call_args)
                    darkimages.append(local_file)
            return darkimages

        if not_null('~{probe_file}'):
            call_args.append('--probe-set=~{probe_file}')

        if not_null('~{target_panel}'):
            call_args.append('--target-panel=~{target_panel}')

        has_image = not_null('~{image}')
        darkimages = get_darkimages('~{sep=";" darkimage}', '~{darkimagestr}')
        has_cimage = not_null('~{colorizedimage}')

        ntrue = has_image + (len(darkimages) > 0) + has_cimage
        if ntrue == 0:
            print("Please set one of the following arguments: image, darkimage or colorizedimage!", file = sys.stderr)
            sys.exit(1)
        elif ntrue > 1:
            print("Please only set one of the following arguments: image, darkimage or colorizedimage!", file = sys.stderr)
            sys.exit(1)

        if has_image:
            call_args.append('--image=~{image}')
        elif len(darkimages) > 0:
            call_args.extend(['--darkimage=' + x for x in darkimages])
        else:
            call_args.append('--colorizedimage=~{colorizedimage}')

        if ('~{area}' is '') and ('~{slide}' is ''):
            call_args.append('--unknown-slide')
        else:
            if '~{area}' is '':
                print("Please provide an input for the 'area' argument!", file = sys.stderr)
                sys.exit(1)
            if '~{slide}' is '':
                print("Please provide an input for the 'slide' argument!", file = sys.stderr)
                sys.exit(1)
            call_args.extend(['--area=~{area}', '--slide=~{slide}'])
            if not_null('~{slidefile}'):
                call_args.append('--slidefile=~{slidefile}')

        if '~{reorient_images}' is 'true':
            call_args.append('--reorient_images')
        if not_null('~{loupe_alignment}'):
            if not has_image:
                print("image option must be set if loupe_alignment is set!", file = sys.stderr)
                sys.exit(1)
            call_args.append('--loupe_alignment=~{loupe_alignment}')

        if '~{no_bam}' is 'true':
            call_args.append('--no-bam')
        if '~{secondary}' is not 'true':
            call_args.append('--nosecondary')

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend ~{backend} -m results/outs "~{output_directory}/~{sample_id}"
    }

    output {
        String output_count_directory = "~{output_directory}/~{sample_id}"
        String output_metrics_summary = "~{output_directory}/~{sample_id}/metrics_summary.csv"
        String output_web_summary = "~{output_directory}/~{sample_id}/web_summary.html"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/spaceranger:~{spaceranger_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
