//
// A nested subworkflow to call structural variants from nuclear DNA.
//

include { CALL_SV_CANVAS                 } from '../../../subworkflows/local/call_sv_canvas/main'
include { CALL_SV_CNVNATOR          } from '../call_sv_cnvnator'
include { CALL_SV_GERMLINECNVCALLER } from '../call_sv_germlinecnvcaller'
include { CALL_SV_MANTA             } from '../call_sv_manta'
include { CALL_SV_TIDDIT            } from '../call_sv_tiddit'

workflow CALL_SV {

    take:
        ch_bwa_index                          // channel: [mandatory] [ val(meta), path(index)]
        ch_canvas_common_cnvs_bed             // channel: [optional] [ val(meta), path(common_cnvs_bed)]
        ch_canvas_f_ploidy_vcf                // channel: [optional, mandatory for canvas] [ val(meta). path(vcf) ]
        ch_canvas_filter_bed                  // channel: [optional] [ val(meta), path(filter13)]
        ch_canvas_kmer_fasta                  // channel: [optional, mandatory for canvas] [val(meta), path(fasta)]
        ch_canvas_m_ploidy_vcf                // channel: [optional, mandatory for canvas] [ val(meta). path(vcf) ]
        ch_case_info                          // channel: [mandatory] [ val(case_info) ]
        ch_gcnvcaller_model                   // channel: [optional; used by mandatory for GATK's cnvcaller][ path(gcnvcaller_model) ]
        ch_genome_bai                         // channel: [mandatory] [ val(meta), path(bai) ]
        ch_genome_bam                         // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_bam_bai                     // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_dictionary                  // channel: [optional; used by mandatory for GATK's cnvcaller][ val(meta), path(dict) ]
        ch_genome_fai                         // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta                       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genomesizes                        // channel: [optional; mandatory for canvas] [ val(meta), path(xml)]
        ch_manta_regions                      // channel: [mandatory] [ path(bed), path(tbi) ]
        ch_ploidy_model                       // channel: [optional; used by mandatory for GATK's cnvcaller][ path(ploidy_model) ]
        ch_readcount_intervals                // channel: [optional; used by mandatory for GATK's cnvcaller][ path(intervals) ]
        ch_snv_vcf                            // channel: [optional; mandatory for canvas] [ val(meta), path(vcf) ]
        ch_target_bed                         // channel: [optional] [ val(meta), path(vcf) ]
        skip_germlinecnvcaller                // boolean
        val_analysis_type                     // string: "wes", "wgs", or "mito"
        val_canvas_reformat_vcf               // boolean: [optional] [default: true] whether to reformat the canvas vcf output to match the expected nf-core/cnvkit vcf output
        val_cnv_caller                        // string: comma-separated list of CNV callers to run. Available callers: canvas and cnvnator

    main:
        ch_canvas_vcf      = channel.empty()
        ch_canvas_seg      = channel.empty()
        ch_cnvnator_vcf    = channel.empty()
        ch_gcnvcaller_vcf  = channel.empty()
        ch_manta_vcf       = channel.empty()
        ch_tiddit_vcf      = channel.empty()

        // CALL_SV is only invoked for non-mito analysis types (gated at the call site in
        // raredisease.nf, mirroring CALL_SV_MT's val_run_mt gate), so no mito check is needed here.
        ch_manta_vcf = CALL_SV_MANTA (ch_genome_bam, ch_genome_bai, ch_genome_fasta, ch_genome_fai, ch_case_info, ch_manta_regions)
            .filtered_diploid_sv_vcf
            .collect{ _meta, vcf -> vcf }

        if (val_analysis_type.equals("wgs")) {
            ch_tiddit_vcf = CALL_SV_TIDDIT (ch_genome_bam_bai, ch_genome_fai, ch_genome_fasta, ch_bwa_index, ch_case_info)
                .vcf
                .collect{ _meta, vcf -> vcf }

            if (val_cnv_caller.contains("cnvnator")) {
                ch_cnvnator_vcf = CALL_SV_CNVNATOR (ch_genome_bam_bai, ch_genome_fasta, ch_case_info)
                    .vcf
                    .collect{ _meta, vcf -> vcf }
            }

            if (val_cnv_caller.contains("canvas")) {
                CALL_SV_CANVAS(
                    ch_genome_bam_bai,
                    ch_snv_vcf,
                    ch_canvas_kmer_fasta,
                    ch_genomesizes,
                    ch_canvas_m_ploidy_vcf,
                    ch_canvas_f_ploidy_vcf,
                    ch_canvas_filter_bed,
                    ch_target_bed,
                    ch_genome_fai,
                    ch_canvas_common_cnvs_bed,
                    val_canvas_reformat_vcf
                )
                ch_canvas_vcf = CALL_SV_CANVAS.out.vcf
                    .collect { _meta, vcf -> vcf }
                ch_canvas_seg = CALL_SV_CANVAS.out.seg
            }

        }

        if (!skip_germlinecnvcaller) {
            ch_gcnvcaller_vcf = CALL_SV_GERMLINECNVCALLER (ch_genome_bam_bai, ch_genome_fasta, ch_genome_fai, ch_readcount_intervals, ch_genome_dictionary, ch_ploidy_model, ch_gcnvcaller_model, ch_case_info)
                .genotyped_filtered_segments_vcf
                .collect{ _meta, vcf -> vcf }

        }

        // Collect individual caller VCFs in a fixed, consistent order so the flat list lines
        // up positionally with the svcaller_priority tags built in main.nf: tiddit -> manta ->
        // gcnvcaller -> cnvnator -> canvas.
        // Empty channels won't contribute any items.
        ch_vcf_paths = ch_tiddit_vcf
            .concat(ch_manta_vcf)
            .concat(ch_gcnvcaller_vcf)
            .concat(ch_cnvnator_vcf)
            .concat(ch_canvas_vcf)
            .collect()

    emit:
        vcfs       = ch_vcf_paths  // channel: [ [path(vcf), path(vcf), ...] ] - flat list of nuclear caller VCFs, in priority order
        canvas_seg = ch_canvas_seg // channel: [ val(meta), path(seg) ]
}
