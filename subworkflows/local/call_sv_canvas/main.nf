include { BEDTOOLS_COMPLEMENT     } from "../../../modules/nf-core/bedtools/complement/main"
include { BEDTOOLS_SUBTRACT       } from "../../../modules/nf-core/bedtools/subtract/main"
include { BCFTOOLS_NORM           } from "../../../modules/nf-core/bcftools/norm/main"
include { BCFTOOLS_VIEW           } from "../../../modules/nf-core/bcftools/view/main"
include { BCFTOOLS_REHEADER       } from "../../../modules/nf-core/bcftools/reheader/main"
include { CANVAS_GERMLINE         } from '../../../modules/local/canvas/germline/main'
include { GAWK as GAWK_CREATE_SEG } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_REFORMAT   } from '../../../modules/nf-core/gawk/main'
include { TABIX_BGZIP             } from '../../../modules/nf-core/tabix/bgzip/main'

workflow CALL_SV_CANVAS {

    take:
    ch_bam_bai // channel: [meta, path(bam), path(bai) ]
    ch_snv     // channel: [meta, path(vcf)]
    ch_kmer_fasta   // channel: [meta, path(fasta)]
    ch_genomesizes // channel: [meta, path(xml)]
    ch_m_ploidy_vcf // channel: [meta, path(male_ploidy_vcf)]
    ch_f_ploidy_vcf // channel: [meta, path(female_ploidy_vcf)]
    ch_canvas_filter_bed // channel: [optional] [meta, path(filter13)]
    ch_target_bed // channel: [optional] [meta, path(filter13)]
    ch_fai // channel: [meta, path(fai)]
    ch_common_cnvs_bed // channel: [optional] [meta, path(common_cnvs_bed)]
    val_reformat_vcf // boolean: [optional] [default: false] whether to reformat the canvas vcf output to match the expected nf-core/cnvkit vcf output

    main:

    // target BED contains regions to include, while Canvas filter contains regions to exclude.
    // So we subtract the filter from the target to get the final regions to include in Canvas.
    ch_subtract_in = ch_target_bed
        .filter { _meta, bed, _index -> bed } // If no target bed is provided, the bed value will be [] (falsey)
        .combine(ch_canvas_filter_bed)
        .map { target_meta, target_bed, _target_index, _filter_meta, canvas_filter_bed ->
            tuple(target_meta, target_bed, canvas_filter_bed)
        }

    // ch_subtract_in will be empty if no target bed is provided, so the subtraction will be skipped implicitily
    BEDTOOLS_SUBTRACT(
       ch_subtract_in
    )

    // Canvas expects a BED file of regions to exclude, so we take the complement of the BED file produced by the subtraction above.
    BEDTOOLS_COMPLEMENT(
        BEDTOOLS_SUBTRACT.out.bed,
        ch_fai.map { _meta, fai -> fai }
    )

    // If no target is provided, the output of BEDTOOLS_COMPLEMENT will be empty, so we use the original canvas filter bed instead
    // In that case, the provided canvas filter bed is the first element in the channel after concatenation
    // Otherwise, the first element is the output of BEDTOOLS_COMPLEMENT
    ch_canvas_filter = BEDTOOLS_COMPLEMENT.out.bed
        .concat(ch_canvas_filter_bed)
        .first()

    // Select correct ploidy VCF based on sex
    ch_ploidy_vcf = ch_bam_bai
        .combine(ch_m_ploidy_vcf)
        .combine(ch_f_ploidy_vcf)
        .map { meta, _bam, _bai, _meta_m, male_ploidy_vcf, _meta_f, female_ploidy_vcf ->
            def ploidy_vcf = meta.sex == 1 ? male_ploidy_vcf : female_ploidy_vcf
            tuple(meta, ploidy_vcf, [], [])
        }

    // Replace placeholder sample names in ploidy VCF with actual sample name
    BCFTOOLS_REHEADER(
        ch_ploidy_vcf,
        [[], []]
    )

    CANVAS_GERMLINE(
        ch_bam_bai,
        ch_kmer_fasta,
        ch_genomesizes,
        ch_canvas_filter,
        ch_snv,
        [[],[]],
        BCFTOOLS_REHEADER.out.vcf,
        ch_common_cnvs_bed
    )

    GAWK_CREATE_SEG(
        CANVAS_GERMLINE.out.covandvarfreq,
        channel.value(file(moduleDir + '/bin/create_seg.awk')),
        false
    )

    ch_bcftools_in = CANVAS_GERMLINE.out.vcf
        .map { meta, vcf -> tuple(meta, vcf, [])}


    if (val_reformat_vcf) {
        BCFTOOLS_NORM(
            ch_bcftools_in,
            ch_kmer_fasta
        )

        GAWK_REFORMAT(
            BCFTOOLS_NORM.out.vcf,
            channel.value(file(moduleDir + '/bin/reformat_canvas.awk')),
            false
        )

        ch_vcf = GAWK_REFORMAT.out.output
            .map { meta, vcf -> tuple(meta, vcf, [])}

    } else {
        ch_vcf = ch_bcftools_in
    }

    BCFTOOLS_VIEW(
        ch_vcf,
        [],
        [],
        []
    )

    emit:
    vcf = BCFTOOLS_VIEW.out.vcf
    seg = GAWK_CREATE_SEG.out.output

}
