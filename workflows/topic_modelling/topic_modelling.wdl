version 1.0

workflow topic_modelling {
  input {
    Int preemptible = 2
    String prefix_exclude = 'mt-,Rpl,Rps'
    Int extra_disk_space = 10
    Float? min_percent_expressed
    Float? max_percent_expressed
    File input_file
    Array[Int] number_of_topics
    Int random_number_seed = 0
    Float disk_space_multiplier = 2
    Int num_cpu = 4
    String memory = "6G"
    String docker = "cumulusprod/topic_modelling:3.8.3"
  }

    call topic_modelling_prepare {
          input:
                preemptible = preemptible,
                prefix_exclude = prefix_exclude,
                disk_space_multiplier = disk_space_multiplier,
                extra_disk_space = extra_disk_space,
                input_file=input_file,
                memory = memory,
                docker = docker
    }

  scatter(topic in number_of_topics) {
    call topic_modelling_gensim {
        input:
              preemptible = preemptible,
              disk_space_multiplier = disk_space_multiplier,
              extra_disk_space = extra_disk_space,
              corpus=topic_modelling_prepare.corpus,
              dictionary=topic_modelling_prepare.dictionary,
              cell_ids=topic_modelling_prepare.cell_ids,
              number_of_topics = topic,
              random_number_seed = random_number_seed,
              num_cpu = num_cpu,
              memory = memory,
              docker = docker
      }
  }
  if(length(number_of_topics)>1) {
      call plot_stats {
        input:
             preemptible = preemptible,
             docker = docker,
             stats = topic_modelling_gensim.stats
      }
  }

  output {
        File? coherence_plot = plot_stats.coherence_plot
        File? perplexity_plot = plot_stats.perplexity_plot
        File corpus = topic_modelling_prepare.corpus
        File dictionary = topic_modelling_prepare.dictionary
        Array[File] stats = topic_modelling_gensim.stats
        Array[File] model = topic_modelling_gensim.model
        Array[File] cell_scores = topic_modelling_gensim.cell_scores
        Array[File] feature_topics = topic_modelling_gensim.feature_topics
        Array[File] report = topic_modelling_gensim.report
  }
}

task plot_stats {
  input {
        Int preemptible
        Array[File] stats
        String docker
  }


  output {
    File monitoring = "monitoring.log"
    File coherence_plot = "coherence.png"
    File perplexity_plot = "perplexity.png"
  }

  command <<<

    set -e

    /software/monitor_script.sh > monitoring.log &

    python /software/lda.py plot \
    ~{sep=' ' stats}
  >>>

  runtime {
    preemptible: preemptible
    bootDiskSizeGb: 12
    disks: "local-disk " + ceil(size(stats, "GB")  + 1) + " HDD"
    docker: "~{docker}"
    cpu: 1
    memory: "1G"
  }

}

task topic_modelling_prepare {
  input {
        Int preemptible
        Int extra_disk_space
        Float disk_space_multiplier
        File input_file
        Float? min_percent_expressed
        Float? max_percent_expressed
        String prefix_exclude
        String memory
        String docker
  }


  output {
    File monitoring = "monitoring.log"
    File corpus = "corpus.mm"
    File dictionary = "dictionary.dict"
    File cell_ids = 'cell_ids.csv'
  }

  command <<<

    set -e

    /software/monitor_script.sh > monitoring.log &

    python /software/lda.py prepare \
    --input ~{input_file} \
    ~{"--prefix_exclude " + prefix_exclude} \
    ~{"--min_percent " + min_percent_expressed} \
    ~{"--max_percent " + max_percent_expressed}
  >>>

  runtime {
    preemptible: preemptible
    bootDiskSizeGb: 12
    disks: "local-disk " + ceil(size(input_file, "GB") * disk_space_multiplier + extra_disk_space) + " HDD"
    docker: "~{docker}"
    cpu: 1
    memory: memory
  }

}

task topic_modelling_gensim {
  input {
        Int preemptible
        Int extra_disk_space
        Float disk_space_multiplier
        File corpus
        File dictionary
        File cell_ids
        Int number_of_topics
        Int random_number_seed
        Int num_cpu
        String memory
        String docker
  }


  output {
    File monitoring = "monitoring.log"
    File model = "lda.model"
    File cell_scores = 'cell_scores.csv'
    File feature_topics = 'feature_topics.csv'
    File stats = 'stats.txt'
    File report = 'report.html'
  }

  command <<<

    set -e

    /software/monitor_script.sh > monitoring.log &

    python /software/lda.py run \
    --cell_ids ~{cell_ids} \
    --corpus ~{corpus} \
    --dictionary ~{dictionary} \
    --topics ~{number_of_topics} \
    --random_seed ~{random_number_seed}
  >>>

  runtime {
    preemptible: preemptible
    bootDiskSizeGb: 12
    disks: "local-disk " + ceil(size(corpus, "GB") * disk_space_multiplier + extra_disk_space) + " HDD"
    docker: "~{docker}"
    cpu: num_cpu
    memory: memory
  }

}

