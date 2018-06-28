workflow cellranger_count {
    # A tab separated file with 2 columns. The first column contains the sample name.
    # The second column contains a the directory containing R1, R2, and I1 fastq files
    # In the following example, sample s1 was sequenced on 2 flow cells while sample s2 was sequenced on 1 flow cell
    # s1 gs://my-bucket/flowcell_1_s1/
    # s1 gs://my-bucket/flowcell_2_s1/
    # s2 gs://my-bucket/flowcell_1_s2/
    File sample_sheet

    # gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz
       # gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz
       # mm10, GRCh38, or a URL to a tar.gz file
       String transcriptome

    File transcriptome_file = (if transcriptome == 'GRCh38'
                                   then 'gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz'
                                   else (if transcriptome == 'mm10'
                                           then 'gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz'
                                           else (if transcriptome == 'GRCh38_and_mm10'
                                                   then 'gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38_and_mm10-1.2.0.tar.gz'
                                                   else transcriptome)))

    # Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization)
    Boolean? secondary
    # Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells
    Int? force_cells
    # Expected number of recovered cells. Default: 3,000 cells. Mutually exclusive with force_cells
    Int? expect_cells
    # Disk space needed per task.
    String disk_space = "250"
    # 2.1.1 or 2.0.2
    String cell_ranger_version = "2.1.1"
    # Number of cpus per cellranger job
    Int cpu = 64
    # Number of preemptible tries
    Int preemptible = 2


    call parse_cellranger_count_input {
        input:
         input_file=sample_sheet,
         preemptible=preemptible,
         cell_ranger_version=cell_ranger_version
    }
    Array[Array[String]] params = parse_cellranger_count_input.params

    scatter (row in params) {
        call CellRangerCount {
            input:
            sample_id = row[0],
            fastq_directory = row[1],
            transcriptome_file = transcriptome_file,
            secondary = secondary,
            expect_cells = expect_cells,
            disk_space = disk_space,
            cell_ranger_version=cell_ranger_version,
            force_cells=force_cells,
            directory_input=true,
            cpu=cpu,
            preemptible=preemptible
       }

   }
}

task CellRangerCount {
    String sample_id
    File transcriptome_file
    String disk_space
    String cell_ranger_version
    Boolean directory_input

    String? fastq_directory
    Int? expect_cells
    Boolean? secondary
    Int? force_cells

    Int cpu
    Int preemptible

    command {
        set -e
        export TMPDIR=/tmp
        python <<CODE
        import os
        from subprocess import call
        transcriptome = "${transcriptome_file}"
        if transcriptome.endswith('.tgz') or transcriptome.endswith('tar.gz'):
            os.mkdir('transcriptome_dir')
            call(['tar', 'xf', transcriptome, '-C', 'transcriptome_dir', '--strip-components', '1'])
        else:
            call(['ln', '-s', transcriptome, 'transcriptome_dir'])
        dirs = dict()
        os.mkdir('fastqs')
        bucket_dirs = '${fastq_directory}'.split(',')
        for i in range(len(bucket_dirs)):
            bucket_dir = bucket_dirs[i]
            os.mkdir('fastqs/' +  str(i))
            call(['gsutil', '-q', '-m', 'cp', '-r', bucket_dir, 'fastqs/' + str(i) + '/'])
        for root, subdirs, files in os.walk('fastqs'):
            has_fastqs = False
            for f in files:
                if f.find('.fastq') > 0:
                    has_fastqs = True
                    break
            if has_fastqs:
                dirs.setdefault(os.path.abspath(root), True)

        keys = list(dirs.keys())
        if len(keys) == 0:
            print('No fastqs found')
            exit(1)
        expect_cells = '${expect_cells}'
        secondary = '${secondary}'
        force_cells = '${force_cells}'
        call_args = list()
        call_args.append('cellranger')
        call_args.append('count')
        call_args.append('--jobmode=local')
        call_args.append('--transcriptome=transcriptome_dir')
        call_args.append('--id=' + '${sample_id}')
        call_args.append('--fastqs=' + ','.join(keys))
        if secondary is not 'true':
            call_args.append('--nosecondary')
        if expect_cells is not '' and force_cells is '':
            call_args.append('--expect-cells=' + str(expect_cells))
        if force_cells is not '':
            call_args.append('--force-cells=' + str(force_cells))
        call(call_args)
        CODE
        }

    output {
        File filtered_gene_bc_matrices_h5 = "${sample_id}/outs/filtered_gene_bc_matrices_h5.h5"
        File metrics_summary = "${sample_id}/outs/metrics_summary.csv"
        File web_summary = "${sample_id}/outs/web_summary.html"
        File possorted_genome_bam = "${sample_id}/outs/possorted_genome_bam.bam"
        File raw_gene_bc_matrices_h5 = "${sample_id}/outs/raw_gene_bc_matrices_h5.h5"

#        File barcodes = "${sample_id}/outs/filtered_gene_bc_matrices/${reference}/barcodes.tsv"
#        File genes = "${sample_id}/outs/filtered_gene_bc_matrices/${reference}/genes.tsv"
#        File matrix = "${sample_id}/outs/filtered_gene_bc_matrices/${reference}/matrix.mtx"

#        File molecule_info = "${sample_id}/outs/molecule_info.h5"
#        File raw_gene_bc_matrices_h5 = "${sample_id}/outs/raw_gene_bc_matrices_h5.h5"

    }

    runtime {
        docker: "regevlab/cellranger-${cell_ranger_version}"
        memory: "416 GB"
        bootDiskSizeGb: 12
        disks: "local-disk ${disk_space} HDD"
        cpu: "${cpu}"
        preemptible: "${preemptible}"
    }
}

task parse_cellranger_count_input {
    File input_file
    Int preemptible
    String cell_ranger_version

    command {
        set -e
        python <<CODE
        import pandas as pd
        df = pd.read_csv('${input_file}', header=None, sep='\t', names=['sample', 'path'])
        df['sample'] = df['sample'].astype(str).str.replace(' ', '')
        samples = df['sample'].unique()
        with open('sample_sheet.txt', 'w') as w:
            for sample in samples:
                paths = df.loc[df['sample'] == sample]['path'].values
                if len(paths) == 0:
                    raise ValueError('No URL for ' + sample)
                w.write(str(sample) + '\t' + ','.join(paths))
                w.write('\n')
        CODE
    }

    output {
        Array[Array[String]] params = read_tsv("sample_sheet.txt")
    }

    runtime {
        preemptible: "${preemptible}"
        docker: "regevlab/cellranger-${cell_ranger_version}"
    }
}
