workflow dropseq {
    File inputSamplesFile
    # tab delimitted text file with name, read1, read2
    Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

    #Reference inputs
    # ["pure"] OR ["HUMAN","MOUSE"]
    Array[String] reference_type = ["pure"]
    File reference_refflat
    File reference_exon_interval
    File reference_gene_interval
    File reference_rrna_interval
    File genome_fasta
    File genome_dict
    File star_genome_refs_zipped

    String five_adaptor = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
    String primer = "AAGCAGTGGTATCAACGCAGAGTAC"
    String? java_memory
    Int collapse_barcodes_threads = 2
    Int star_threads = 6
    Int base_qc = 10
    Int min_core_reads = 3000
    Int core_barcodes = 1000
    # use core_barcodes value directly instead of estimating from elbow plot
    Boolean force_cells = false
    Int min_noncore_reads = 30000000
    # maximum available RAM for sorting BAM
    String limitBAMsortRAM = "45000000000"
    # runtime node memory
    String star_memory = "100 GB"
    # max number of collapsed junctions
    Int limitOutSJcollapsed = 1000000
    Int max_records = 10000000
    Int edit_distance = 1
    Int umi_edit_distance = 1
    Int molecular_edit_distance = 1
    Int cb_low_quality_base_fail = 1
    Int umi_bases_fail_quality = 1
    Int collapse_read_quality = 10
    Int read_map_quality = 10
    Int polya_mismatch = 0
    Int polya_min_base = 6
    Int smart_mismatches = 0
    Int smart_min_base = 5
    Int max_synthesis_errors = 2
    Int plot_point_size = 1

    String collapsed_cell_tag = "XC"
    String cell_tag = "XC"
    String gene_exon_tag = "GE"
    String molecular_tag = "XM"
    String polya_trim_tag = "ZP"
    String start_trim_tag = "ZS"

    String intermediate_dir = "intermediate_dir"
    String temp = "temp"

  scatter (sample in inputSamples) {
          call FastqToSam {
              input:
                   java_memory=java_memory,
                   fastq1=sample[1],
                   fastq2=sample[2],
                   max_records=max_records,
                   base_qc=base_qc,
                   sample_name=sample[0],
                   temp=temp,
                   intermediate_dir=intermediate_dir
          }

          call TagBamWithReadSequenceExtended as TagXM {
              input:
                   java_memory=java_memory,
                   input_sam=FastqToSam.raw_Sam,
                   tagged_sam=sample[0],
                   start_base=13,
                   stop_base=20,
                   base_qc=base_qc,
                   tag=molecular_tag,
                   umi_bases_fail_quality=umi_bases_fail_quality,
                   discard_read="False",
                   summary="unaligned_tagged1.bam_summary",
                   temp=temp
                                                       }


          call TagBamWithReadSequenceExtended as TagXC {
              input:
                   java_memory=java_memory,
                   input_sam=TagXM.Tagged_Sam,
                   tagged_sam=sample[0],
                   start_base=1,
                   stop_base=12,
                   base_qc=base_qc,
                   tag=cell_tag,
                   umi_bases_fail_quality=cb_low_quality_base_fail,
                   discard_read="True",
                   summary="unaligned_tagged3.bam_summary",
                   temp=temp
                                                       }

          call FilterXQBam {
              input:

                  java_memory=java_memory,
                  xc_tagged_sam=TagXC.Tagged_Sam,
                  tag=cell_tag,
                  base_qc=base_qc,
                  sample_name=sample[0],
                  temp=temp,
                  intermediate_dir=intermediate_dir
          }

          call TrimStartingSequence {
              input:

                  java_memory=java_memory,
                  tagged_filtered_sam=FilterXQBam.Tagged_Filtered_Sam,
                  five_adaptor=five_adaptor,
                  smart_mismatches=smart_mismatches,
                  smart_min_base=smart_min_base,
                  base_qc=base_qc,
                  sample_name=sample[0],
                  temp=temp,
                  intermediate_dir=intermediate_dir
                                    }

           call PolyATrimmer {
               input:

                   java_memory=java_memory,
                   trimmed_sam=TrimStartingSequence.Trimmed_Sam,
                   polya_mismatch=polya_mismatch,
                   polya_min_base=polya_min_base,
                   polya_trim_summary="polyA_trimming_report.txt",
                   polya_trim_unaligned_bam="unaligned_mc_tagged_polyA_filtered.bam",
                   temp=temp
                              }

           call SamToFastq {
               input:

                   java_memory=java_memory,
                   polya_trim_unaligned_bam=PolyATrimmer.Polya_Trim_Unaligned_Bam,
                   polya_trim_unaligned_fq="unaligned_tagged.fastq.gz",
                   temp=temp
            }

           call star_align {
               input:
                    star_genome_refs_zipped=star_genome_refs_zipped,
                    star_threads=star_threads,
                    align_fastq=SamToFastq.Polya_Trim_Unaligned_Fq,
                    limitBAMsortRAM=limitBAMsortRAM,
                    star_memory=star_memory,
                    sample_name=sample[0],
                    base_name="star",
                    limitOutSJcollapsed=limitOutSJcollapsed,
                    base_qc=base_qc,
                    temp=temp
            }

            call MergeBamAlignment {
                input:
                     java_memory=java_memory,
                     genome_dict=genome_dict,
                     genome_fasta=genome_fasta,
                     polya_trim_unaligned_bam=PolyATrimmer.Polya_Trim_Unaligned_Bam,
                     star_bam=star_align.Output_Bam,
                     max_records=max_records,
                     intermediate_dir=intermediate_dir
                                   }


            call TagReadWithInterval as TagXG {
                input:
                     java_memory=java_memory,
                     in_bam=MergeBamAlignment.Merged_Bam,
                     feat_intervals=reference_gene_interval,
                     tag="XG",
                     intermediate_dir=intermediate_dir
                                              }


            call TagReadWithInterval as TagXE {
                input:

                     java_memory=java_memory,
                     in_bam=TagXG.Tag_Interval_Out,
                     feat_intervals=reference_exon_interval,
                     tag="XE",
                     intermediate_dir=intermediate_dir
                                             }

            call TagReadWithGeneExon {
                input:

                      java_memory=java_memory,
                      tag_interval_exon=TagXE.Tag_Interval_Out,
                      ref_flat=reference_refflat,
                      exon_tag=gene_exon_tag,
                      base_qc=base_qc,
                      sample_name=sample[0],
                      temp=temp
                                     }

            call CollapseBarCodes {
                input:

                 java_memory=java_memory,
                 tag_read=TagReadWithGeneExon.Tag_Read,
                 cell_tag=cell_tag,
                 edit_distance=edit_distance,
                 collapse_read_quality=collapse_read_quality,
                 max_records=max_records,
                 collapse_barcodes_threads=collapse_barcodes_threads,
                 min_core_read=min_core_reads,
                 min_noncore_read=min_noncore_reads,
                 base_qc=base_qc,
                 sample_name=sample[0]
                                  }

            call DetectBeadSynthesisErrors {
                input:

                 java_memory=java_memory,
                 collapsed_barcodes=CollapseBarCodes.Collapsed_Barcodes,
                 core_barcode=2*core_barcodes,
                 primer=primer,
                 max_synthesis_errors=max_synthesis_errors,
                 base_qc=base_qc,
                 sample_name=sample[0]
                                           }

            call SamIndex {
                input:
                 corrected_barcodes=DetectBeadSynthesisErrors.Corrected_Barcodes
                          }

            call BAMTagHistogram {
                input:

                 java_memory=java_memory,
                 corrected_barcodes=DetectBeadSynthesisErrors.Corrected_Barcodes,
                 cell_tag_collapsed=collapsed_cell_tag,
                 base_qc=base_qc,
                 sample_name=sample[0]
                                 }

            call decompress_hist {
                input:
                 bam_tag_histogram=BAMTagHistogram.Bam_Tag_Histogram,
                 base_qc=base_qc,
                 sample_name=sample[0]
                                 }

            call DropSeqCumuPlot {
                input:
                 bam_tag_histogram_unzip=decompress_hist.Bam_Tag_Histogram_Unzip,
                 sample_name=sample[0],
                 base_qc=base_qc
                                 }

            call CollectCellBarcodes {
                input:
                 bam_tag_histogram_unzip=decompress_hist.Bam_Tag_Histogram_Unzip,
                 drop_seq_cum_num_cells=if force_cells then core_barcodes else if force_cells then core_barcodes else read_int(DropSeqCumuPlot.Drop_Seq_Cum_Num_Cells_Text),
                 sample_name=sample[0],
                 base_qc=base_qc
                                     }

           call CollectRnaSeqMetrics {
               input:
                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  reflat=reference_refflat,
                  ribo_interval=reference_rrna_interval,
                  java_memory=java_memory,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                     }

           call SingleCellRnaSeqMetricsCollector {
               input:
                  java_memory=java_memory,
                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  reflat=reference_refflat,
                  ribo_interval=reference_rrna_interval,
                  cell_tag_collapsed=collapsed_cell_tag,
                  estimated_cells=if force_cells then core_barcodes else read_int(DropSeqCumuPlot.Drop_Seq_Cum_Num_Cells_Text),
                  read_map_quality=read_map_quality,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                                 }

           call GatherReadQualityMetrics {
               input:
                  java_memory=java_memory,

                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                         }

           call GatherReadQualityMetricsByCell {
               input:
                  java_memory=java_memory,
                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  cell_tag_collapsed=collapsed_cell_tag,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                               }

           call MeanQualityByCycleAllReads {
               input:
                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  java_memory=java_memory,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                           }

           call MeanQualityPerCycleAlignedReads {
               input:
                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  java_memory=java_memory,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                                }

           call TrimTag as StartTrimTag {
               input:
                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  java_memory=java_memory,
                  tag=start_trim_tag,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                        }

           call TrimTag as PolyATrimTag {
               input:
                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  java_memory=java_memory,
                  tag=polya_trim_tag,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                        }

           call NumReadsPerCellBarcode as NumReadsPerCellBarcodeCell {
               input:
                  java_memory=java_memory,
                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  tag=collapsed_cell_tag,
                  reads_per_cell_barcode_quality=10,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam"),
                  output_filename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")+"_"+cell_tag+"_TrimHist.txt"
                                                                      }


           call NumReadsPerCellBarcode as NumReadsPerCellBarcodeCellCollapsed {
               input:
                  java_memory=java_memory,
                  input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                  tag=collapsed_cell_tag,
                  reads_per_cell_barcode_quality=10,
                  file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam"),
                  output_filename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")+"_numReads_perCell_"+cell_tag+"_mq_10_collapsed.txt.gz"
                                                                              }

            call BarcodeBaseDistribution as BarcodeBaseDistributionCellTag {
                input:
                   java_memory=java_memory,
                   input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                   tag=cell_tag,
                   file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                                                           }

            call BarcodeBaseDistribution as BarcodeBaseDistributionMolTag {
                input:
                   java_memory=java_memory,
                   input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                   tag=molecular_tag,
                   file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                                                          }
            call selectCellsByReadsSCM {
                input:
                   input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                   reads_per_cell_barcode_col=NumReadsPerCellBarcodeCellCollapsed.Num_Reads_Per_Cell_Barcode,
                   estimated_cells=if force_cells then core_barcodes else read_int(DropSeqCumuPlot.Drop_Seq_Cum_Num_Cells_Text),
                   file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                       }

            call copyBeadSynthesisErrorDetail {
                input:
                   bead_error=DetectBeadSynthesisErrors.Bead_Error_Out,
                   bead_error_basename=basename(DetectBeadSynthesisErrors.Bead_Error_Out)
                                              }

            if (length(reference_type) == 2) {

                  call FilterMixedBams as FilterMixedBamsDigitalExpression {
                      input:

                         java_memory=java_memory,
                         corrected_barcodes=DetectBeadSynthesisErrors.Corrected_Barcodes,
                         reference_type=reference_type,
                         sample_name=sample[0],
                         out_dir="bams"
                                                                           }

                  call MixedDigitalExpressionUMI {
                      input:
                         java_memory=java_memory,
                         filter_inputs=FilterMixedBamsDigitalExpression.Output_filtered,
                         umi_edit_distance=umi_edit_distance,
                         cell_tag_collapsed=collapsed_cell_tag,
                         barcodes=CollectCellBarcodes.Barcodes,
                         sample_name=sample[0]
                                                 }
                   call MixedDigitalExpressionReads {
                       input:
                          java_memory=java_memory,
                          filter_inputs=FilterMixedBamsDigitalExpression.Output_filtered,
                          umi_edit_distance=umi_edit_distance,
                          cell_tag_collapsed=collapsed_cell_tag,
                          barcodes=CollectCellBarcodes.Barcodes,
                          sample_name=sample[0]
                                                     }

                   call FilterMixedBams as FilterMixedBamsforSummarize {
                       input:
                        java_memory=java_memory,
                        corrected_barcodes=DetectBeadSynthesisErrors.Corrected_Barcodes,
                        reference_type=reference_type,
                        sample_name=sample[0],
                        out_dir="reports"
                                                                   }

                   call DigitalExpressionSpeciesMultiple {
                       input:
                        java_memory=java_memory,
                        input_bams=FilterMixedBamsforSummarize.Output_filtered,
                        file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam"),
                        reference_type=reference_type,
                        collapsed_cell_tag=collapsed_cell_tag,
                        molecular_tag=molecular_tag,
                        gene_exon_tag=gene_exon_tag,
                        molecular_edit_distance=molecular_edit_distance,
                        read_map_quality=read_map_quality,
                        selected_barcodes=selectCellsByReadsSCM.Selected_Barcodes
                                                       }

                  call MolecularBarcodeDistributionsByGeneMultiple {
                      input:
                       java_memory=java_memory,
                       input_bams=FilterMixedBamsforSummarize.Output_filtered,
                       file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam"),
                       collapsed_cell_tag=collapsed_cell_tag,
                       molecular_tag=molecular_tag,
                       gene_exon_tag=gene_exon_tag,
                       molecular_edit_distance=molecular_edit_distance,
                       read_map_quality=read_map_quality,
                       selected_barcodes=selectCellsByReadsSCM.Selected_Barcodes
                                                                 }

                 call categorizeCellsUsingKnee {
                     input:
                       file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam"),
                       knee_docs=DigitalExpressionSpeciesMultiple.Species_Digital_Summary,
                       selected_barcodes=selectCellsByReadsSCM.Selected_Barcodes,
                       species_list=reference_type,
                       estimated_cells=if force_cells then core_barcodes else read_int(DropSeqCumuPlot.Drop_Seq_Cum_Num_Cells_Text)
                                               }

                 call NumReadsPerCellBarcodeOrganism {
                     input:
                       java_memory=java_memory,
                       species_bams=FilterMixedBamsforSummarize.Output_filtered,
                       cell_tag_collapsed=collapsed_cell_tag,
                       read_map_quality=read_map_quality,
                       sample_name=sample[0]
                                                     }

                 call plotPairOrganism {
                     input:
                     qc_metrics=GatherReadQualityMetrics.QC_metrics,
                     cycle_quality_metrics=MeanQualityByCycleAllReads.Cycle_Quality_Metrics,
                     cycle_quality_metrics_aligned=MeanQualityPerCycleAlignedReads.Cycle_Quality_Metrics_Aligned,
                     frac_intron_exon_file=CollectRnaSeqMetrics.Frac_Intron_Exon_File,
                     frac_intron_exon_cell_file=SingleCellRnaSeqMetricsCollector.Frac_Intron_Exon_Cell_File,
                     start_hist=StartTrimTag.TrimTagHist,
                     polya_hist=PolyATrimTag.TrimTagHist,
                     reads_per_cell_barcode=NumReadsPerCellBarcodeCell.Num_Reads_Per_Cell_Barcode,
                     reads_per_cell_barcode_col=NumReadsPerCellBarcodeCellCollapsed.Num_Reads_Per_Cell_Barcode,
                     estimated_cells=if force_cells then core_barcodes else read_int(DropSeqCumuPlot.Drop_Seq_Cum_Num_Cells_Text),
                     selected_barcodes=selectCellsByReadsSCM.Selected_Barcodes,
                     barcode_distribution_cell=BarcodeBaseDistributionCellTag.Barcode_Distribution,
                     barcode_distribution_mol=BarcodeBaseDistributionMolTag.Barcode_Distribution,
                     cell_qc_metrics=GatherReadQualityMetricsByCell.Cell_QC_Metrics,
                     point_size=plot_point_size,
                     species_digital_summary=DigitalExpressionSpeciesMultiple.Species_Digital_Summary,
                     species_barcode_file=MolecularBarcodeDistributionsByGeneMultiple.Species_Barcode_File,
                     updated_bead_error=copyBeadSynthesisErrorDetail.Bead_Error_Copy,
                     file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam"),
                     Species_Reads_Cell=NumReadsPerCellBarcodeOrganism.Species_Reads_Cell,
                     knee_output=categorizeCellsUsingKnee.knee_output
                                       }
                                            }

            if (length(reference_type) == 1 ) {
                call PureDigitalExpressionUMI {
                      input:
                         java_memory=java_memory,
                         filter_input=DetectBeadSynthesisErrors.Corrected_Barcodes,
                         umi_edit_distance=umi_edit_distance,
                         cell_tag_collapsed=collapsed_cell_tag,
                         barcodes=CollectCellBarcodes.Barcodes,
                         sample_name=sample[0],
                         reference_type=reference_type[0]
                                               }
                call PureDigitalExpressionReads {
                      input:

                         java_memory=java_memory,
                         filter_input=DetectBeadSynthesisErrors.Corrected_Barcodes,
                         umi_edit_distance=umi_edit_distance,
                         cell_tag_collapsed=collapsed_cell_tag,
                         barcodes=CollectCellBarcodes.Barcodes,
                         sample_name=sample[0],
                         reference_type=reference_type[0]
                                                 }

                call DigitalExpressionSpeciesSingle {
                    input:
                      java_memory=java_memory,
                      input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                      file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam"),
                      reference_type=reference_type[0],
                      collapsed_cell_tag=collapsed_cell_tag,
                      molecular_tag=molecular_tag,
                      gene_exon_tag=gene_exon_tag,
                      molecular_edit_distance=molecular_edit_distance,
                      read_map_quality=read_map_quality,
                      selected_barcodes=selectCellsByReadsSCM.Selected_Barcodes
                }

                call MolecularBarcodeDistributionsByGeneSingle {
                    input:
                      java_memory=java_memory,
                      input_bam=DetectBeadSynthesisErrors.Corrected_Barcodes,
                      file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam"),
                      collapsed_cell_tag=collapsed_cell_tag,
                      molecular_tag=molecular_tag,
                      gene_exon_tag=gene_exon_tag,
                      molecular_edit_distance=molecular_edit_distance,
                      read_map_quality=read_map_quality,
                      selected_barcodes=selectCellsByReadsSCM.Selected_Barcodes
                }

                call plotSingleOrganism {
                    input:
                     qc_metrics=GatherReadQualityMetrics.QC_metrics,
                     cycle_quality_metrics=MeanQualityByCycleAllReads.Cycle_Quality_Metrics,
                     cycle_quality_metrics_aligned=MeanQualityPerCycleAlignedReads.Cycle_Quality_Metrics_Aligned,
                     frac_intron_exon_file=CollectRnaSeqMetrics.Frac_Intron_Exon_File,
                     frac_intron_exon_cell_file=SingleCellRnaSeqMetricsCollector.Frac_Intron_Exon_Cell_File,
                     start_hist=StartTrimTag.TrimTagHist,
                     polya_hist=PolyATrimTag.TrimTagHist,
                     reads_per_cell_barcode=NumReadsPerCellBarcodeCell.Num_Reads_Per_Cell_Barcode,
                     reads_per_cell_barcode_col=NumReadsPerCellBarcodeCellCollapsed.Num_Reads_Per_Cell_Barcode,
                     estimated_cells=if force_cells then core_barcodes else read_int(DropSeqCumuPlot.Drop_Seq_Cum_Num_Cells_Text),
                     selected_barcodes=selectCellsByReadsSCM.Selected_Barcodes,
                     barcode_distribution_cell=BarcodeBaseDistributionCellTag.Barcode_Distribution,
                     barcode_distribution_mol=BarcodeBaseDistributionMolTag.Barcode_Distribution,
                     cell_qc_metrics=GatherReadQualityMetricsByCell.Cell_QC_Metrics,
                     point_size=plot_point_size,
                     species_digital_single_summary=DigitalExpressionSpeciesSingle.Species_Digital_Summary,
                     species_barcode_single_file=MolecularBarcodeDistributionsByGeneSingle.Species_Barcode_File,
                     updated_bead_error=copyBeadSynthesisErrorDetail.Bead_Error_Copy,
                     file_basename=basename(DetectBeadSynthesisErrors.Corrected_Barcodes,".bam")
                                         }
                }
     }
}

task FastqToSam {
 String? java_memory
 File fastq1
 File fastq2
 String intermediate_dir
 Int max_records
 String sample_name
 String temp
 Int base_qc

 command {
   java -Xmx${default="4g" java_memory} \
   -jar /home/picard.jar \
   FastqToSam \
   FASTQ=${fastq1} \
   FASTQ2=${fastq2} \
   TMP_DIR=${temp}/${intermediate_dir}_${base_qc}_${sample_name} \
   QUALITY_FORMAT=Standard \
   OUTPUT=${temp}/${intermediate_dir}_${base_qc}_${sample_name}/${sample_name}.sam \
   COMPRESSION_LEVEL=0 \
   MAX_RECORDS_IN_RAM=${max_records} \
   SAMPLE_NAME=${sample_name} \
   SORT_ORDER=queryname
 }

 output {
    File raw_Sam="${temp}/${intermediate_dir}_${base_qc}_${sample_name}/${sample_name}.sam"
 }

 runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub((((size(fastq1,"GB")+size(fastq2,"GB"))+1)*20),"\\..*","") + " HDD"
    memory: "50 GB"
    preemptible: 2
 }
}

task TagBamWithReadSequenceExtended {

 String? java_memory
 File input_sam
 String tagged_sam
 String summary
 Int start_base
 Int stop_base
 Int base_qc
 String tag
 Int umi_bases_fail_quality
 String discard_read
 String temp

 command {
   mkdir temp

   TagBamWithReadSequenceExtended \
   -m ${default="4g" java_memory} \
   INPUT=${input_sam} \
   OUTPUT=${temp}/${tagged_sam}_${tag}.sam \
   COMPRESSION_LEVEL=0 \
   SUMMARY=${temp}/${summary}_${tag}.txt \
   BASE_RANGE=${start_base}-${stop_base} \
   BASE_QUALITY=${base_qc} \
   BARCODED_READ=1 \
   DISCARD_READ=${discard_read} \
   TAG_NAME=${tag} \
   NUM_BASES_BELOW_QUALITY=${umi_bases_fail_quality}
 }

 output {
   File Tagged_Sam="${temp}/${tagged_sam}_${tag}.sam"
   File Summary="${temp}/${summary}_${tag}.txt"
 }

  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_sam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task FilterXQBam {

 String? java_memory
 File xc_tagged_sam
 String tag
 String temp
 String intermediate_dir
 String base_qc
 String sample_name

 command {
   mkdir temp
   mkdir ${temp}/${intermediate_dir}_${base_qc}_${sample_name}

   FilterBAM \
   -m ${default="4g" java_memory} \
   TAG_REJECT=XQ \
   INPUT=${xc_tagged_sam} \
   OUTPUT=${temp}/${intermediate_dir}_${base_qc}_${sample_name}/${sample_name}_${tag}_filtered.sam \
   COMPRESSION_LEVEL=0
 }

 output {
   File Tagged_Filtered_Sam="${temp}/${intermediate_dir}_${base_qc}_${sample_name}/${sample_name}_${tag}_filtered.sam"
 }

  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(xc_tagged_sam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task TrimStartingSequence {

 String? java_memory
 File tagged_filtered_sam
 String five_adaptor
 Int smart_mismatches
 Int smart_min_base
 String intermediate_dir
 String base_qc
 String sample_name
 String temp

 command {
   mkdir temp
   mkdir temp/${intermediate_dir}_${base_qc}_${sample_name}

   TrimStartingSequence \
   -m ${default="4g" java_memory} \
   INPUT=${tagged_filtered_sam} \
   OUTPUT=${temp}/${intermediate_dir}_${base_qc}_${sample_name}/${sample_name}_filtered_trim.sam \
   OUTPUT_SUMMARY=${temp}/adaptor_trimming_report.txt \
   SEQUENCE=${five_adaptor} \
   MISMATCHES=${smart_mismatches} \
   NUM_BASES=${smart_min_base}
 }

 output {
  File Trimmed_Sam="${temp}/${intermediate_dir}_${base_qc}_${sample_name}/${sample_name}_filtered_trim.sam"
  File Trim_Summary="${temp}/adaptor_trimming_report.txt"
 }

  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(tagged_filtered_sam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task PolyATrimmer {

 String? java_memory
 File trimmed_sam
 String polya_trim_unaligned_bam
 String polya_trim_summary
 Int polya_mismatch
 Int polya_min_base
 String temp

 command {
   mkdir temp

   PolyATrimmer \
   -m ${default="4g" java_memory} \
   INPUT=${trimmed_sam} \
   OUTPUT=${temp}/${polya_trim_unaligned_bam} \
   OUTPUT_SUMMARY=${temp}/${polya_trim_summary} \
   MISMATCHES=${polya_mismatch} \
   NUM_BASES=${polya_min_base}
 }

 output {
  File Polya_Trim_Unaligned_Bam="${temp}/${polya_trim_unaligned_bam}"
  File Polya_Trim_Summary="${temp}/${polya_trim_summary}"
 }

  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(trimmed_sam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task SamToFastq {

 String? java_memory
 File polya_trim_unaligned_bam
 String polya_trim_unaligned_fq
 String temp

 command {
   mkdir temp

   java -Xmx${default="4g" java_memory} \
   -jar /home/picard.jar \
   SamToFastq \
   INPUT=${polya_trim_unaligned_bam} \
   FASTQ=${temp}/${polya_trim_unaligned_fq} \
   tmp_dir=${temp} \
   CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2
 }

 output {
  File Polya_Trim_Unaligned_Fq="${temp}/${polya_trim_unaligned_fq}"
 }

  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(polya_trim_unaligned_bam,"GB")+1)*4),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task star_align {

 Int star_threads
 File align_fastq
 String limitBAMsortRAM
 String base_name
 String sample_name
 String temp
 String base_qc
 String limitOutSJcollapsed
 String star_memory
 File star_genome_refs_zipped

 command {

 tar xvzf ${star_genome_refs_zipped}

  STAR --genomeDir STAR2_5 \
  --runThreadN ${star_threads} \
  --readFilesIn ${align_fastq} \
  --readFilesCommand "gunzip -c" \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \
  --limitBAMsortRAM ${limitBAMsortRAM} \
  --outFileNamePrefix ${base_name} \
  --limitOutSJcollapsed ${limitOutSJcollapsed}
 }

 output {
  File Output_Bam="${base_name}Aligned.sortedByCoord.out.bam"
  File Output_Log_Final="${base_name}Log.final.out"
  File Output_Log="${base_name}Log.out"
  File Output_Log_Progress="${base_name}Log.progress.out"
  File Output_SJ="${base_name}SJ.out.tab"
 }

   runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(align_fastq,"GB")+30)*20),"\\..*","") + " HDD"
    memory: "${star_memory}"
    preemptible: 2
    cpu: "${star_threads}"
 }

}

task MergeBamAlignment {
 String intermediate_dir
 String? java_memory

 File genome_fasta
 File genome_dict
 File polya_trim_unaligned_bam
 File star_bam
 Int max_records

 command {
  mkdir ${intermediate_dir}

  java -Djava.io.tmpdir=${intermediate_dir} \
  -Xmx${default="4g" java_memory} -jar \
  /home/picard.jar MergeBamAlignment \
  REFERENCE_SEQUENCE=${genome_fasta} \
  UNMAPPED_BAM=${polya_trim_unaligned_bam} \
  ALIGNED_BAM=${star_bam} \
  MAX_RECORDS_IN_RAM=${max_records} \
  OUTPUT=${intermediate_dir}/starAligned.sortedByCoord.out_merged.bam \
  COMPRESSION_LEVEL=0 \
  INCLUDE_SECONDARY_ALIGNMENTS=false \
  PAIRED_RUN=False VALIDATION_STRINGENCY=SILENT
 }

 output {
  File Merged_Bam="${intermediate_dir}/starAligned.sortedByCoord.out_merged.bam"
 }

 runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(polya_trim_unaligned_bam,"GB")+30)*20),"\\..*","") + " HDD"
    memory: "100GB"
    preemptible: 2
 }
}

task TagReadWithInterval {

 String? java_memory
 File in_bam
 File feat_intervals
 String tag
 String intermediate_dir

 command {
  mkdir ${intermediate_dir}

  TagReadWithInterval \
  -m ${default="4g" java_memory} \
  I=${in_bam} \
  O=${intermediate_dir}/starAligned.sortedByCoord.out_merged_${tag}.bam \
  COMPRESSION_LEVEL=0 \
  LOCI=${feat_intervals} \
  TAG=${tag} || : ;
 }

 output {
  File Tag_Interval_Out="${intermediate_dir}/starAligned.sortedByCoord.out_merged_${tag}.bam"
 }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(in_bam,"GB")+30)*20),"\\..*","") + " HDD"
    memory: "100GB"
    preemptible: 2
 }
}

task TagReadWithGeneExon {

 String? java_memory
 File tag_interval_exon
 File ref_flat
 String exon_tag
 String sample_name
 String base_qc
 String temp

 command {
  mkdir -p ${temp}/star_fast_${base_qc}_${sample_name}

  TagReadWithGeneExon \
  -m ${default="4g" java_memory} \
  I=${tag_interval_exon} \
  O=${temp}/star_fast_${base_qc}_${sample_name}/${sample_name}_star_gene_exon_tagged2.bam \
  ANNOTATIONS_FILE=${ref_flat} \
  TAG=${exon_tag}
 }

 output {
  File Tag_Read="${temp}/star_fast_${base_qc}_${sample_name}/${sample_name}_star_gene_exon_tagged2.bam"
 }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(tag_interval_exon,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task CollapseBarCodes {

 String? java_memory
 File tag_read
 String cell_tag
 Int edit_distance
 Int collapse_read_quality
 Int max_records
 Int collapse_barcodes_threads
 Int min_core_read
 Int min_noncore_read
 String sample_name
 Int base_qc

 command {
  mkdir "bams"

  CollapseBarcodesInPlace \
  -m ${default="4g" java_memory} \
  INPUT=${tag_read} \
  OUTPUT=bams/${sample_name}_bq${base_qc}_star.bam \
  PRIMARY_BARCODE=${cell_tag} \
  OUT_BARCODE=ZC \
  EDIT_DISTANCE=${edit_distance} \
  READ_QUALITY=${collapse_read_quality} \
  FILTER_PCR_DUPLICATES=false \
  FIND_INDELS=true \
  MAX_RECORDS_IN_RAM=${max_records} \
  NUM_THREADS=${collapse_barcodes_threads} \
  MIN_NUM_READS_CORE=${min_core_read} \
  MIN_NUM_READS_NONCORE=${min_noncore_read}
 }

 output {
  File Collapsed_Barcodes="bams/${sample_name}_bq${base_qc}_star.bam"
 }
   runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(tag_read,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    cpu: "${collapse_barcodes_threads}"
    preemptible: 2
 }
}

task DetectBeadSynthesisErrors {

 String? java_memory
 File collapsed_barcodes
 Int core_barcode
 String primer
 Int max_synthesis_errors
 String sample_name
 Int base_qc

 command {
  mkdir "bams"
  mkdir "bam_reads"

  DetectBeadSynthesisErrors \
    -m ${default="4g" java_memory} \
    I=${collapsed_barcodes} \
    O=bams/${sample_name}_bq${base_qc}_star_corrected.bam \
    OUTPUT_STATS=bam_reads/${sample_name}_bq${base_qc}.synthesis_stats.bead_synthesis_error_detail \
    SUMMARY=bam_reads/${sample_name}_bq${base_qc}.synthesis_stats.summary.txt \
    NUM_BARCODES=${core_barcode} \
    PRIMER_SEQUENCE=${primer} \
    MAX_NUM_ERRORS=${max_synthesis_errors}
 }

 output {
  File Corrected_Barcodes="bams/${sample_name}_bq${base_qc}_star_corrected.bam"
  File Bead_Error_Out="bam_reads/${sample_name}_bq${base_qc}.synthesis_stats.bead_synthesis_error_detail"
  File Bead_Error_Summary="bam_reads/${sample_name}_bq${base_qc}.synthesis_stats.summary.txt"
 }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(collapsed_barcodes,"GB")+5)*2),"\\..*","") + " HDD"
    memory: "100GB"
    preemptible: 2
 }

}

task SamIndex {
 File corrected_barcodes

 command {
  samtools index \
  ${corrected_barcodes}
 }

 output {
  File Corrected_Barcodes="${corrected_barcodes}.bai"
 }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(corrected_barcodes,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task BAMTagHistogram {

 String? java_memory
 File corrected_barcodes
 String cell_tag_collapsed
 Int base_qc
 String sample_name

 command {
  mkdir "bam_reads"

  BAMTagHistogram \
  -m ${default="4g" java_memory} \
  I=${corrected_barcodes} \
  O=bam_reads/${sample_name}_bq${base_qc}_star.reads.txt.gz \
  TAG=${cell_tag_collapsed}
 }

 output {
  File Bam_Tag_Histogram="bam_reads/${sample_name}_bq${base_qc}_star.reads.txt.gz"
 }
   runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(corrected_barcodes,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task decompress_hist {
 File bam_tag_histogram
 Int base_qc
 String sample_name

 command {
  mkdir "bam_reads"
  gunzip -c ${bam_tag_histogram} > \
  bam_reads/${sample_name}_bq${base_qc}_star.reads.txt
 }

 output {
  File Bam_Tag_Histogram_Unzip="bam_reads/${sample_name}_bq${base_qc}_star.reads.txt"
 }
   runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(bam_tag_histogram,"GB")+1)*2),"\\..*","") + " HDD"
    preemptible: 2
 }
}

task DropSeqCumuPlot {
 File bam_tag_histogram_unzip
 String sample_name
 Int base_qc

 command {
  mkdir "bam_reads"

  DropSeqCumuPlot.R \
  --collapsed ${bam_tag_histogram_unzip} \
  --counts ${bam_tag_histogram_unzip} \
  --num_cells bam_reads/${sample_name}_bq${base_qc}_star_numCells.txt \
  --cummulative_plot bam_reads/${sample_name}_bq${base_qc}_star_cumplot.pdf \
  --reads_plot bam_reads/${sample_name}_bq${base_qc}_star_NreadsHiToLo.pdf
 }

 output {
  File Drop_Seq_Cum_Num_Cells_Text="bam_reads/${sample_name}_bq${base_qc}_star_numCells.txt"
  File Drop_Seq_Cum_Plots="bam_reads/${sample_name}_bq${base_qc}_star_cumplot.pdf"
  File Drop_Seq_Read_Plots="bam_reads/${sample_name}_bq${base_qc}_star_NreadsHiToLo.pdf"
 }
 runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(bam_tag_histogram_unzip,"GB")+1)*2),"\\..*","") + " HDD"
    preemptible: 2
 }
}

task CollectCellBarcodes {
 File bam_tag_histogram_unzip
 Int drop_seq_cum_num_cells

 String base_qc
 String sample_name

 command {
  mkdir "bam_reads"

  collect_cell_barcodes.R \
  ${bam_tag_histogram_unzip} \
  ${drop_seq_cum_num_cells} \
  bam_reads/${sample_name}_bq${base_qc}_star_star.barcodes_use.txt
 }

 output {
  File Barcodes="bam_reads/${sample_name}_bq${base_qc}_star_star.barcodes_use.txt"
 }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(bam_tag_histogram_unzip,"GB")+1)*2),"\\..*","") + " HDD"
    preemptible: 2
 }
}

task FilterMixedBams {

 String? java_memory
 File corrected_barcodes
 Array[String] reference_type
 String sample_name
 String out_dir

 command {
  mkdir ${out_dir}

  for key in ${sep=" " reference_type}
  do
      FilterBAM \
      -m ${default="4g" java_memory} \
      I=${corrected_barcodes} \
      O=${out_dir}/${sample_name}.$key.bam \
      REF_SOFT_MATCHED_RETAINED=$key
  done
 }

 output {
  Array[Pair[File,String]] Output_filtered=[("${out_dir}/${sample_name}.${reference_type[0]}.bam", "${reference_type[0]}"), ("${out_dir}/${sample_name}.${reference_type[1]}.bam", "${reference_type[1]}")]
 }
 runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(corrected_barcodes,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task PureDigitalExpressionUMI {

 String? java_memory
 File filter_input
 Int umi_edit_distance
 String cell_tag_collapsed
 File barcodes
 String sample_name
 String reference_type

 command {
  mkdir "UMI_DGE"

    DigitalExpression \
    -m ${default="4g" java_memory} \
    I=${filter_input} \
    OUTPUT=UMI_DGE/${sample_name}_${reference_type}.umi.dge.txt.gz \
    EDIT_DISTANCE=${umi_edit_distance} \
    CELL_BARCODE_TAG=${cell_tag_collapsed} \
    SUMMARY=UMI_DGE/${sample_name}_${reference_type}.umi.dge.summary.txt \
    CELL_BC_FILE=${barcodes}
         }

  output {
    File Dge_Umi_Output="UMI_DGE/${sample_name}_${reference_type}.umi.dge.txt.gz"
    File Dge_Umi_Summary="UMI_DGE/${sample_name}_${reference_type}.umi.dge.summary.txt"
  }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(filter_input,"GB")+size(barcodes,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}


task MixedDigitalExpressionUMI {

 String? java_memory
 Array[Pair[File,String]] filter_inputs
 Int umi_edit_distance
 String cell_tag_collapsed
 File barcodes
 String sample_name

 command {
  mkdir "UMI_DGE"

  python <<CODE
  import os
  from subprocess import call
  for file_name, species in [${sep="," filter_inputs}]:
      call(["DigitalExpression",
            "-m","${default="4g" java_memory}",
            "I="+file_name,
            "OUTPUT="+os.path.join("UMI_DGE","${sample_name}"+"_"+species+".umi.dge.txt.gz"),
            "EDIT_DISTANCE="+str(${umi_edit_distance}),
            "CELL_BARCODE_TAG="+"${cell_tag_collapsed}",
            "SUMMARY="+os.path.join("UMI_DGE","${sample_name}"+"_"+species+".umi.dge.summary.txt"),
            "CELL_BC_FILE="+"${barcodes}"])
  CODE
  }

  output {
    Array[File] Dge_Umi_Output=glob("UMI_DGE/${sample_name}_*.umi.dge.txt.gz")
    Array[File] Dge_Umi_Summary=glob("UMI_DGE/${sample_name}_*.umi.dge.summary.txt")
  }

  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(filter_inputs[0].left,"GB")+size(barcodes,"GB")+1)*3),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}

task PureDigitalExpressionReads {

 String? java_memory
 File filter_input
 Int umi_edit_distance
 String cell_tag_collapsed
 File barcodes
 String sample_name
 String reference_type

 command {
  mkdir "reads_DGE"

    DigitalExpression \
    -m ${default="4g" java_memory} \
    I=${filter_input} \
    OUTPUT=reads_DGE/${sample_name}_${reference_type}.reads.dge.txt.gz  \
    EDIT_DISTANCE=${umi_edit_distance} \
    CELL_BARCODE_TAG=${cell_tag_collapsed} \
    SUMMARY=reads_DGE/${sample_name}_${reference_type}.reads.dge.summary.txt \
    CELL_BC_FILE=${barcodes} \
    OUTPUT_READS_INSTEAD=true
  }

  output {
      File Dge_Reads_Output="reads_DGE/${sample_name}_${reference_type}.reads.dge.txt.gz"
      File Dge_Reads_Summary="reads_DGE/${sample_name}_${reference_type}.reads.dge.summary.txt"
  }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(filter_input,"GB")+size(barcodes,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
 }
}


task MixedDigitalExpressionReads {

 String? java_memory
 Array[Pair[File,String]] filter_inputs
 Int umi_edit_distance
 String cell_tag_collapsed
 File barcodes
 String sample_name

 command {
  mkdir "reads_DGE"

  python <<CODE
  import os
  from subprocess import call
  for file_name, species in [${sep="," filter_inputs}]:
      call(["DigitalExpression",
            "-m","${default="4g" java_memory}",
            "I="+file_name,
            "OUTPUT="+os.path.join("reads_DGE","${sample_name}"+"_"+species+".reads.dge.txt.gz"),
            "EDIT_DISTANCE="+str(${umi_edit_distance}),
            "CELL_BARCODE_TAG="+"${cell_tag_collapsed}",
            "SUMMARY="+os.path.join("reads_DGE","${sample_name}"+"_"+species+".reads.dge.summary.txt"),
            "CELL_BC_FILE="+"${barcodes}",
            "OUTPUT_READS_INSTEAD=true"])
  CODE
  }

  output {
      Array[File] Dge_Reads_Output=glob("reads_DGE/${sample_name}_*.reads.dge.txt.gz")
      Array[File] Dge_Reads_Summary=glob("reads_DGE/${sample_name}_*.reads.dge.summary.txt")
  }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(filter_inputs[0].left,"GB")+size(barcodes,"GB")+1)*3),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
}

 task CollectRnaSeqMetrics {
  File input_bam
  File reflat
  File ribo_interval
  String? java_memory

  String file_basename

  command {
   mkdir "reports"

   java -Xmx${default="4g" java_memory} \
   -jar /home/picard.jar CollectRnaSeqMetrics \
   I=${input_bam} \
   REF_FLAT=${reflat} \
   STRAND_SPECIFICITY=NONE \
   OUTPUT=reports/${file_basename}_fracIntronicExonic.txt \
   RIBOSOMAL_INTERVALS=${ribo_interval}
          }

  output {
   File Frac_Intron_Exon_File="reports/${file_basename}_fracIntronicExonic.txt"
         }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+2)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
}

  task SingleCellRnaSeqMetricsCollector {
   String? java_memory
   File input_bam
   File reflat
   File ribo_interval
   String cell_tag_collapsed
   Int estimated_cells
   Int read_map_quality

   String file_basename

  command {
    mkdir "reports"

    SingleCellRnaSeqMetricsCollector \
    -m ${default="4g" java_memory} \
    I=${input_bam} \
    ANNOTATIONS_FILE=${reflat} \
    OUTPUT=reports/${file_basename}_fracIntronicExonicPerCell.txt.gz \
    RIBOSOMAL_INTERVALS=${ribo_interval} \
    CELL_BARCODE_TAG=${cell_tag_collapsed} \
    NUM_CORE_BARCODES=${estimated_cells} \
    READ_MQ=${read_map_quality}
  }

  output {
    File Frac_Intron_Exon_Cell_File="reports/${file_basename}_fracIntronicExonicPerCell.txt.gz"
  }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+2)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
}


task GatherReadQualityMetrics {
   String? java_memory
   File input_bam

   String file_basename

  command {
    mkdir "reports"

    GatherReadQualityMetrics \
    -m ${default="4g" java_memory} \
    I=${input_bam} \
    O=reports/${file_basename}_ReadQualityMetrics.txt
  }

  output {
    File QC_metrics="reports/${file_basename}_ReadQualityMetrics.txt"
  }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
}

task GatherReadQualityMetricsByCell {
 String? java_memory
 File input_bam
 String cell_tag_collapsed

 String file_basename

 command {
  mkdir "reports"

  GatherReadQualityMetrics \
  -m ${default="4g" java_memory} \
  I=${input_bam} \
  O=reports/${file_basename}_ReadQualityMetricsByCell.txt.gz \
  TAG=${cell_tag_collapsed}
         }

 output {
   File Cell_QC_Metrics="reports/${file_basename}_ReadQualityMetricsByCell.txt.gz"
        }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
                                    }

 task MeanQualityByCycleAllReads {
  File input_bam
  String? java_memory

  String file_basename

  command {
   mkdir "reports"

   java -Xmx${default="4g" java_memory} \
   -jar /home/picard.jar \
   MeanQualityByCycle \
   I=${input_bam} \
   OUTPUT=reports/${file_basename}_meanQualityPerCycle_allReads.txt \
   CHART=reports/${file_basename}_meanQualityPerCycle_allReads.pdf
          }

  output {
   File Cycle_Quality_Metrics="reports/${file_basename}_meanQualityPerCycle_allReads.txt"
   File Cycle_Quality_Metrics_PDF="reports/${file_basename}_meanQualityPerCycle_allReads.pdf"
         }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }

                                }

 task MeanQualityPerCycleAlignedReads {
  File input_bam
  String? java_memory

  String file_basename

  command {
   mkdir "reports"

   java -Xmx${default="4g" java_memory} \
   -jar /home/picard.jar \
   MeanQualityByCycle \
   I=${input_bam} \
   OUTPUT=reports/${file_basename}_meanQualityPerCycle_alignedReads.txt \
   CHART=reports/${file_basename}_meanQualityPerCycle_alignedReads.pdf \
   ALIGNED_READS_ONLY=true
          }

  output {
   File Cycle_Quality_Metrics_Aligned="reports/${file_basename}_meanQualityPerCycle_alignedReads.txt"
   File Cycle_Quality_Metrics_Aligned_PDF="reports/${file_basename}_meanQualityPerCycle_alignedReads.pdf"
         }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
                                      }

 task TrimTag {
  File input_bam
  String? java_memory
  String tag

  String file_basename

  command {
   mkdir "reports"

   BAMTagHistogram \
   -m ${default="4g" java_memory} \
   I=${input_bam} \
   OUTPUT=reports/${file_basename}_${tag}_TrimHist.txt \
   TAG=${tag} \
   FILTER_PCR_DUPLICATES=false \
   READ_QUALITY=0
          }
  output {
     File TrimTagHist="reports/${file_basename}_${tag}_TrimHist.txt"
         }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
              }

 task NumReadsPerCellBarcode {
  String? java_memory
  File input_bam
  String tag
  Int reads_per_cell_barcode_quality
  String file_basename
  String output_filename

  command {
   mkdir "reports"

   BAMTagHistogram \
   -m ${default="4g" java_memory} \
   I=${input_bam} \
   OUTPUT=reports/${output_filename} \
   TAG=${tag} \
   FILTER_PCR_DUPLICATES=false \
   READ_QUALITY=${reads_per_cell_barcode_quality}
          }

  output {
   File Num_Reads_Per_Cell_Barcode="reports/${output_filename}"
         }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
                            }

  task BarcodeBaseDistribution {
   String? java_memory
   File input_bam
   String tag

   String file_basename

   command {
    mkdir "reports"

    BaseDistributionAtReadPosition \
    -m ${default="4g" java_memory} \
    I=${input_bam} \
    OUTPUT=reports/${file_basename}_barcode_distribution_${tag}.txt \
    TAG=${tag}
           }

   output {
     File Barcode_Distribution="reports/${file_basename}_barcode_distribution_${tag}.txt"
          }
    runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
                                }

   task selectCellsByReadsSCM {
    File input_bam
    File reads_per_cell_barcode_col
    String file_basename
    Int estimated_cells
    Int? estimated_beads
    Int beadcount = if defined(estimated_beads) then select_first([estimated_beads]) else (estimated_cells * 10)

    command {
     mkdir "reports"

     echo $R_LIBS

     Rscript -e 'library(DropSeq.barnyard)' -e 'selectCellsByReadsSCM(bamFile="${input_bam}",reportDir="reports",readsPerCellBarcodeFile="${reads_per_cell_barcode_col}",estimatedNumBeads=${beadcount},outputFile="reports/${file_basename}_auto.selectedCellBarcodes.txt")'
    }

    output {
     File Selected_Barcodes="reports/${file_basename}_auto.selectedCellBarcodes.txt"
           }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+size(reads_per_cell_barcode_col,"GB")+1)*2),"\\..*","") + " HDD"
    preemptible: 2
  }
}

 task copyBeadSynthesisErrorDetail {
  File bead_error
  String bead_error_basename

  command {
   mkdir "reports"
   cp ${bead_error} reports/${bead_error_basename}
  }

  output {
   File Bead_Error_Copy="reports/${bead_error_basename}"
  }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(bead_error,"GB")+1)*2),"\\..*","") + " HDD"
    preemptible: 2
  }
 }

 task filterBamByOrganism {
  String? java_memory
  File input_bam
  Array[String] reference_type
  String file_basename


  command {
   mkdir "reports"

   for species in ${sep=" " reference_type}
   do
     FilterBAM \
     -m ${default="4g" java_memory} \
     I=${input_bam} \
     O=reports/${file_basename}.$species.bam \
     REF_SOFT_MATCHED_RETAINED=$species
   done
  }

  output {
   Array[Pair[File, String]] Species_Bam = [("reports/${file_basename}.${reference_type[0]}.bam", "${reference_type[0]}"), ("reports/${file_basename}.${reference_type[1]}.bam", "${reference_type[1]}")]
  }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
 }

 task DigitalExpressionSpeciesSingle {
  String? java_memory
  File input_bam
  String file_basename
  String reference_type
  String collapsed_cell_tag
  String molecular_tag
  String gene_exon_tag
  Int molecular_edit_distance
  Int read_map_quality
  File selected_barcodes


  command {
      mkdir "reports"

      DigitalExpression \
        -m ${default="4g" java_memory} \
        I=${input_bam} \
        O=reports/${file_basename}_${reference_type}_auto_digital_expression.txt.gz \
        SUMMARY=reports/${file_basename}_${reference_type}_auto_digital_expression_summary.txt \
        CELL_BARCODE_TAG=${collapsed_cell_tag} \
        MOLECULAR_BARCODE_TAG=${molecular_tag} \
        GENE_EXON_TAG=${gene_exon_tag} \
        EDIT_DISTANCE=${molecular_edit_distance} \
        READ_MQ=${read_map_quality} \
        MIN_BC_READ_THRESHOLD=0 \
        CELL_BC_FILE=${selected_barcodes}
           }

  output {
   File Species_Digital_Expression="reports/${file_basename}_${reference_type}_auto_digital_expression.txt.gz"
   File Species_Digital_Summary="reports/${file_basename}_${reference_type}_auto_digital_expression_summary.txt"
  }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
 }


 task DigitalExpressionSpeciesMultiple {
  String? java_memory
  Array[Pair[File,String]] input_bams
  String file_basename
  Array[String] reference_type
  String collapsed_cell_tag
  String molecular_tag
  String gene_exon_tag
  Int molecular_edit_distance
  Int read_map_quality
  File selected_barcodes


  command {
      mkdir "reports"

      python <<CODE
      import os
      from subprocess import call
      for species_file,species_name in [${sep="," input_bams}]:
          call(["DigitalExpression",
                "-m", "${default="4g" java_memory}",
                "I="+species_file,
                "O="+os.path.join("reports","${file_basename}"+"_"+species_name+"_auto_digital_expression.txt.gz"),
                "SUMMARY="+os.path.join("reports","${file_basename}"+"_"+species_name+"_auto_digital_expression_summary.txt"),
                "CELL_BARCODE_TAG="+"${collapsed_cell_tag}",
                "MOLECULAR_BARCODE_TAG="+"${molecular_tag}",
                "GENE_EXON_TAG="+"${gene_exon_tag}",
                "EDIT_DISTANCE="+str(${molecular_edit_distance}),
                "READ_MQ="+str(${read_map_quality}),
                "MIN_BC_READ_THRESHOLD=0",
                "CELL_BC_FILE="+"${selected_barcodes}"])
      CODE
           }

  output {
   Array[File] Species_Digital_Expression=glob("reports/${file_basename}_*_auto_digital_expression.txt.gz")
   Array[File] Species_Digital_Summary=glob("reports/${file_basename}_*_auto_digital_expression_summary.txt")
  }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bams[0].left,"GB")+1)*4),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
 }

task MolecularBarcodeDistributionsByGeneSingle {

  String? java_memory
  File input_bam
  String file_basename
  String collapsed_cell_tag
  String molecular_tag
  String gene_exon_tag
  Int molecular_edit_distance
  Int read_map_quality
  File selected_barcodes


  command {
     mkdir "reports"

     GatherMolecularBarcodeDistributionByGene \
        -m ${default="4g" java_memory} \
        I=${input_bam} \
        O=reports/${file_basename}_auto_molBC.txt.gz \
        CELL_BARCODE_TAG=${collapsed_cell_tag} \
        MOLECULAR_BARCODE_TAG=${molecular_tag} \
        GENE_EXON_TAG=${gene_exon_tag} \
        EDIT_DISTANCE=${molecular_edit_distance} \
        READ_MQ=${read_map_quality} \
        MIN_BC_READ_THRESHOLD=0 \
        CELL_BC_FILE=${selected_barcodes}
          }

  output {
   File Species_Barcode_File="reports/${file_basename}_auto_molBC.txt.gz"
         }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
                                               }


task MolecularBarcodeDistributionsByGeneMultiple {

  String? java_memory
  Array[Pair[File,String]] input_bams
  String file_basename
  String collapsed_cell_tag
  String molecular_tag
  String gene_exon_tag
  Int molecular_edit_distance
  Int read_map_quality
  File selected_barcodes


  command {
     mkdir "reports"

     python <<CODE
     import os
     from subprocess import call
     for species_bam,species_name in [${sep="," input_bams}]:
         call(["GatherMolecularBarcodeDistributionByGene",
               "-m","${default="4g" java_memory}",
               "I="+species_bam,
               "O="+os.path.join("reports","${file_basename}"+"_"+species_name+"_auto_molBC.txt.gz"),
               "CELL_BARCODE_TAG="+"${collapsed_cell_tag}",
               "MOLECULAR_BARCODE_TAG="+"${molecular_tag}",
               "GENE_EXON_TAG="+"${gene_exon_tag}",
               "EDIT_DISTANCE="+str(${molecular_edit_distance}),
               "READ_MQ="+str(${read_map_quality}),
               "MIN_BC_READ_THRESHOLD=0",
               "CELL_BC_FILE="+"${selected_barcodes}"])
     CODE
          }

  output {
   Array[File] Species_Barcode_File=glob("reports/${file_basename}_*_auto_molBC.txt.gz")
         }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(input_bams[0].left,"GB")+1)*4),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }
                                                  }

task categorizeCellsUsingKnee {

     String file_basename
     Array[File] knee_docs
     Array[String] species_list
     File selected_barcodes

     Int estimated_cells
     Int? estimated_beads
     Int beadcount = if defined(estimated_beads) then select_first([estimated_beads]) else (estimated_cells * 10)

     command {
      mkdir "reports"

      Rscript -e 'library(DropSeq.barnyard)' -e \
        'categorizeCellsUsingKnee(digitalExpressionFile1="${knee_docs[0]}",digitalExpressionFile2="${knee_docs[1]}",organismOne="${species_list[0]}",organismTwo="${species_list[1]}",selectedCellsFile="${selected_barcodes}",outFile="reports/${file_basename}auto_categorized_cellTypes.txt")'
             }

     output {
      File knee_output="reports/${file_basename}auto_categorized_cellTypes.txt"
            }
  runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(selected_barcodes,"GB")+1)*2),"\\..*","") + " HDD"
    preemptible: 2
  }
}

task NumReadsPerCellBarcodeOrganism {
     String? java_memory
     Array[Pair[File,String]] species_bams
     String cell_tag_collapsed
     Int read_map_quality

     String sample_name

     command {
      mkdir "reports"

      python <<CODE
      import os
      from subprocess import call
      for org_bam,species in [${sep="," species_bams}]:
          call(["BAMTagHistogram",
                "-m", "${default="4g" java_memory}",
                "I="+org_bam,
                "OUTPUT="+os.path.join("reports",os.path.splitext(os.path.basename(org_bam))[0]+"_"+"numReads_perCell_"+"${cell_tag_collapsed}"+"_mq_"+str(${read_map_quality})+".txt.gz"),
                "TAG="+"${cell_tag_collapsed}",
                "READ_QUALITY="+str(${read_map_quality}),
                "FILTER_PCR_DUPLICATES=false"])
      CODE
             }

     output {
       Array[Pair[File, String]] Species_Reads_Cell=[("reports/${sample_name}.${species_bams[0].right}_numReads_perCell_${cell_tag_collapsed}_mq_${read_map_quality}.txt.gz","${species_bams[0].right}"),("reports/${sample_name}.${species_bams[1].right}_numReads_perCell_${cell_tag_collapsed}_mq_${read_map_quality}.txt.gz","${species_bams[1].right}")]
            }
      runtime {
    docker: "regevlab/dropseq_v5"
    disks: "local-disk " + sub(((size(species_bams[0].left,"GB")+1)*4),"\\..*","") + " HDD"
    memory: "4GB"
    preemptible: 2
  }

}
task plotSingleOrganism {

     File qc_metrics
     File cycle_quality_metrics
     File cycle_quality_metrics_aligned
     File frac_intron_exon_file
     File frac_intron_exon_cell_file
     File start_hist
     File polya_hist
     File reads_per_cell_barcode
     File reads_per_cell_barcode_col
     Int estimated_cells
     File selected_barcodes
     File barcode_distribution_cell
     File barcode_distribution_mol
     File cell_qc_metrics
     Int point_size
     File species_digital_single_summary
     File species_barcode_single_file
     File updated_bead_error
     String file_basename

     Int? estimated_beads
     Int beadcount = if defined(estimated_beads) then select_first([estimated_beads]) else (estimated_cells * 10)

   command {
     mkdir -p "reports"

     Rscript -e 'library(DropSeq.barnyard)' -e \
     'plotSingleOrganism(outPlot="reports/${file_basename}_auto.pdf",alignmentQualityFile="${qc_metrics}",meanQualityAllFile="${cycle_quality_metrics}",meanQualityAlignedFile="${cycle_quality_metrics_aligned}",exonIntronFile="${frac_intron_exon_file}",exonIntronPerCellFile="${frac_intron_exon_cell_file}",startTagTrimFile="${start_hist}",polyATagTrimFile="${polya_hist}",cellBCCountsFile="${reads_per_cell_barcode}",readsPerCellBarcodeFile="${reads_per_cell_barcode_col}",estimatedNumCells=${estimated_cells},estimatedNumBeads=${beadcount},selectedCellsFile="${selected_barcodes}",basePctMatrixCellFile="${barcode_distribution_cell}",basePctMatrixMolecularFile="${barcode_distribution_mol}",alignmentQualityByCellFile="${cell_qc_metrics}",point.cex=${point_size}, digitalExpressionSummaryFile="${species_digital_single_summary}", molecularBarcodeDistributionByGeneFile="${species_barcode_single_file}",beadSynthesisErrorDetailFile="${updated_bead_error}")'
           }

   output {
      File Pdf="reports/${file_basename}_auto.pdf"
          }
  runtime {
    docker: "regevlab/dropseq_v5"
    preemptible: 2
  }
}

task plotPairOrganism {

     File qc_metrics
     File cycle_quality_metrics
     File cycle_quality_metrics_aligned
     File frac_intron_exon_file
     File frac_intron_exon_cell_file
     File start_hist
     File polya_hist
     File reads_per_cell_barcode
     File reads_per_cell_barcode_col
     Int estimated_cells
     Int? estimated_beads
     Int beadcount = if defined(estimated_beads) then select_first([estimated_beads]) else (estimated_cells * 10)
     File selected_barcodes
     File barcode_distribution_cell
     File barcode_distribution_mol
     File cell_qc_metrics
     Int point_size
     Array[File] species_digital_summary
     Array[File] species_barcode_file
     Array[Pair[File,String]] Species_Reads_Cell
     File updated_bead_error
     String file_basename

     File knee_output

   command {
     mkdir -p "reports/${file_basename}"

     Rscript -e 'library(DropSeq.barnyard)' -e \
     'plotPairOrganism(outPlot="reports/${file_basename}_auto.pdf",alignmentQualityFile="${qc_metrics}",meanQualityAllFile="${cycle_quality_metrics}",meanQualityAlignedFile="${cycle_quality_metrics_aligned}",exonIntronFile="${frac_intron_exon_file}",exonIntronPerCellFile="${frac_intron_exon_cell_file}",startTagTrimFile="${start_hist}",polyATagTrimFile="${polya_hist}",cellBCCountsFile="${reads_per_cell_barcode}",readsPerCellBarcodeFile="${reads_per_cell_barcode_col}",estimatedNumCells=${estimated_cells},estimatedNumBeads=${beadcount},selectedCellsFile="${selected_barcodes}",basePctMatrixCellFile="${barcode_distribution_cell}",basePctMatrixMolecularFile="${barcode_distribution_mol}",alignmentQualityByCellFile="${cell_qc_metrics}",point.cex=${point_size},digitalExpressionSummaryFile1="${species_digital_summary[0]}",digitalExpressionSummaryFile2="${species_digital_summary[1]}",digitalExpressionSummaryFile=NULL,readsPerCellBCOrganismFile1="${Species_Reads_Cell[0].left}",readsPerCellBCOrganismFile2="${Species_Reads_Cell[1].left}",organism1="${Species_Reads_Cell[0].right}",organism2="${Species_Reads_Cell[1].right}",cellTypesFile="${knee_output}",molecularBarcodeDistributionByGeneFile1="${species_barcode_file[0]}",molecularBarcodeDistributionByGeneFile2="${species_barcode_file[1]}",beadSynthesisErrorDetailFile="${updated_bead_error}")'
           }

   output {
      File Pdf="reports/${file_basename}_auto.pdf"
          }
    runtime {
    docker: "regevlab/dropseq_v5"
    preemptible: 2
  }
}


