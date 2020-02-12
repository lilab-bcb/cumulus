version 1.0

import "star-solo.wdl" as sts
import "bustools.wdl" as kbc
import "alevin.wdl" as ale

workflow count {
    input {
        String sample_id
        File r1_fastq
        File r2_fastq
        File? i1_fastq

        String genome
        String chemistry
        String output_directory
        String count_tool

        String docker_registry = "cumulusprod"
        Int num_cpu = 32
        Int disk_space = 100
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        Int memory = 32

        # Star-Solo
        String starsolo_star_version = "2.7.3a"

        # Alevin
        String alevin_version = '1.1'

        # Bustools
        Boolean bustools_output_loom = false
        Boolean bustools_output_h5ad = false
        String bustools_docker = "shaleklab/kallisto-bustools"
        String bustools_version = '0.24.4'
    }
    
    File ref_index_file = "gs://regev-lab/resources/count_tools/ref_index.tsv"
    # File ref_index_file = "ref_index.tsv"
    Map[String, String] ref_index2gsurl = read_map(ref_index_file)
    String genome_url = ref_index2gsurl[genome]

    if (count_tool == 'StarSolo') {
        call sts.starsolo as StarSolo {
            input:
                sample_id = sample_id,
                r1_fastq = r1_fastq,
                r2_fastq = r2_fastq,
                genome_url = genome_url + '/starsolo.tar.gz',
                chemistry = chemistry,
                output_directory = output_directory,
                num_cpu = num_cpu,
                star_version = starsolo_star_version,
                docker_registry = docker_registry,
                disk_space = disk_space,
                preemptible = preemptible,
                zones = zones,
                memory = memory
        }
    }

    if (count_tool == 'Alevin') {
        call ale.alevin as Alevin {
            input:
                sample_id = sample_id,
                r1_fastq = r1_fastq,
                r2_fastq = r2_fastq,
                genome_url = genome_url + '/alevin.tar.gz',
                chemistry = chemistry,
                output_directory = output_directory,
                num_cpu = num_cpu,
                docker_registry = docker_registry,
                alevin_version = alevin_version,
                disk_space = disk_space,
                preemptible = preemptible,
                zones = zones,
                memory = memory
        }
    }

    if (count_tool == 'Bustools') {
        call kbc.bustools as Bustools {
            input:
                sample_id = sample_id,
                r1_fastq = r1_fastq,
                r2_fastq = r2_fastq,
                genome_url = genome_url + '/bustools.tar.gz',
                chemistry = chemistry,
                output_directory = output_directory,
                num_cpu = num_cpu,
                output_loom = bustools_output_loom,
                output_h5ad = bustools_output_h5ad,
                bustools_version = bustools_version,
                disk_space = disk_space,
                preemptible = preemptible,
                zones = zones,
                memory = memory
        }
    }

    output {
        File? starsolo_monitor = StarSolo.monitoringLog
        String? starsolo_output_folder = StarSolo.output_folder
        File? alevin_monitor = Alevin.monitoringLog
        String? alevin_output_folder = Alevin.output_folder
        String? bustools_output_folder = Bustools.output_folder
    }    
}