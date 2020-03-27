workflow demuxlet {
    # 4 column tsv file with no header with name, bam, barcodes, vcf
    File tsv_file
    Array[Array[String]] input_tsv = read_tsv(tsv_file)
    String? zones = "us-central1-b us-east1-d us-west1-a us-west1-b"
    Int? num_cpu = 1
    String? memory = "10G"
    Int? extra_disk_space = 2
    Int? preemptible = 2
    String? docker = "cumulusprod/demuxlet:0.1-beta"
    # gs://fc-xx/out/filtered_feature_bc_matrix.h5
    # gs://fc-xx/out/possorted_genome_bam.bam
    # gs://fc-xx/out/possorted_genome_bam.bam.bai
    # gs://fc-xx/out/filtered_feature_bc_matrix/barcodes.tsv.gz
    Int min_MQ = 20
    Array[Float] alpha = [0.1, 0.2, 0.3, 0.4, 0.5]
    Int min_TD = 0
    String tag_group = "CB"
    String tag_UMI  = "UB"
    String field = "GT"
    Float geno_error = 0.1

    scatter (row in input_tsv) {
        call demuxlet_task {
               input:
                   alpha = alpha,
                   min_TD = min_TD,
                   min_MQ = min_MQ,
                   prefix=row[0],
                   bam = row[1],
                   barcodes=row[2],
                   vcf = row[3],
                   tag_group = tag_group,
                   tag_UMI  = tag_UMI,
                   field = field,
                   geno_error=geno_error,
                   zones=zones,
                   num_cpu=num_cpu,
                   memory=memory,
                   extra_disk_space=extra_disk_space,
                   preemptible=preemptible,
                   docker=docker
           }
    }
}

task demuxlet_task {
    String docker
    File bam
    File vcf
    Array[Float] alpha
    Int min_TD
    Int min_MQ
    Int num_cpu
    String memory
    Int extra_disk_space
    Int preemptible
    String zones
    File barcodes
    String tag_group
    String tag_UMI
    String field
    String prefix
    Float geno_error

     command {
            set -e
            monitor_script.sh > monitoring.log &

            python <<CODE
            import pandas as pd

            df = pd.read_csv('${barcodes}', sep='\t', header=None)
            #if remove_barcode_suffix:
            #     df[0] = df[0].str[:-2]
            df.to_csv('__barcodes.txt', sep='\t', index=False, header=False)

            CODE

            popscle demuxlet --sam ${bam} \
            --tag-group ${tag_group} \
            --tag-UMI ${tag_UMI} \
            --vcf ${vcf} \
            --field ${field} \
            --geno-error-offset ${geno_error} \
            --out ${prefix} \
            --alpha ${sep=' --alpha ' alpha} \
            --min-TD ${min_TD} \
            --min-MQ ${min_MQ} \
            --group-list __barcodes.txt

        }

      output {
        File best  = "${prefix}.best"
        File monitoringLog = "monitoring.log"
    }

    runtime {
           docker: docker
           zones: zones
           memory: memory
           bootDiskSizeGb: 12
           disks: "local-disk " + ceil(extra_disk_space + size(bam,"GB") +  size(vcf,"GB"))  + " HDD"
           cpu: num_cpu
           preemptible: preemptible
       }
}

