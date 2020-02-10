#import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:count_tools_merge_fastqs/versions/2/plain-WDL/descriptor" as mfw
#import "https://raw.githubusercontent.com/HumanCellAtlas/skylab/optimus_v1.4.0_terra/pipelines/optimus/Optimus.wdl" as opm
#import "https://api.firecloud.org/ga4gh/v1/tools/alexandria:kallisto-bustools_count/versions/1/plain-WDL/descriptor" as kbc

import "star-solo.wdl" as sts
#import "alevin.wdl" as ale

workflow count {
    File input_csv_file
    String genome
    String chemistry
    String output_directory
    String count_tool

    String? docker_registry = "cumulusprod"
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? memory = 32

    File ref_index_file = "gs://regev-lab/resources/count_tools/ref_index.tsv"
    # File ref_index_file = "ref_index.tsv"
    Map[String, String] ref_index2gsurl = read_map(ref_index_file)
	String genome_url = ref_index2gsurl[genome]

	# Star-Solo
	Int? starsolo_num_cpu = 32
	String? starsolo_star_version = "2.7.3a"


    if (count_tool == 'StarSolo') {
        call sts.starsolo as star_solo {
            input:
                input_sample_sheet = input_csv_file,
                genome_url = genome_url,
                chemistry = chemistry,
                output_directory = output_directory,
                num_cpu = starsolo_num_cpu,
                star_version = starsolo_star_version,
                docker_registry = docker_registry,
                disk_space = disk_space,
                preemptible = preemptible,
                zones = zones,
                memory = memory
        }
    }

    output {
        File? starsolo_monitor = star_solo.monitoringLog
        String? starsolo_output_folder = star_solo.output_folder
    }    
}