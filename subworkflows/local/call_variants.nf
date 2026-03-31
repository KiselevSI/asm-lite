// subworkflows/local/call_variants.nf

include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM    } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW    } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX   } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_STATS   } from '../../modules/nf-core/bcftools/stats/main'

workflow CALL_VARIANTS {
    take:
    bam_bai_ref  // tuple(meta[+ref_name], bam, bai, ref_fasta, ref_fai)

    main:

    // 1. mpileup + call (combined in nf-core module)
    // BCFTOOLS_MPILEUP(tuple(meta, bam, intervals_mpileup, intervals_call), tuple(meta2, fasta, fai), save_mpileup)
    ch_mpileup_input = bam_bai_ref.map { meta, bam, bai, fasta, fai ->
        [meta, bam, [], []]  // no intervals
    }
    ch_mpileup_ref = bam_bai_ref.map { meta, bam, bai, fasta, fai ->
        [[id: meta.ref_name], fasta, fai]
    }

    BCFTOOLS_MPILEUP(ch_mpileup_input, ch_mpileup_ref, false)

    // 2. norm
    // BCFTOOLS_NORM(tuple(meta, vcf, tbi), tuple(meta2, fasta))
    ch_norm_input = BCFTOOLS_MPILEUP.out.vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }
    ch_norm_ref = bam_bai_ref.map { meta, bam, bai, fasta, fai ->
        [[id: meta.ref_name], fasta]
    }

    BCFTOOLS_NORM(ch_norm_input, ch_norm_ref)

    // 3. view (filter: QUAL>=20 && FMT/DP>=10 via ext.args in modules.config)
    // BCFTOOLS_VIEW(tuple(meta, vcf, index), regions, targets, samples)
    ch_view_input = BCFTOOLS_NORM.out.vcf.map { meta, vcf -> [meta, vcf, []] }
    BCFTOOLS_VIEW(ch_view_input, [], [], [])

    // 4. index
    BCFTOOLS_INDEX(BCFTOOLS_VIEW.out.vcf)

    // 5. stats
    // BCFTOOLS_STATS(tuple(meta, vcf, tbi), tuple(meta2,regions), tuple(meta3,targets), tuple(meta4,samples), tuple(meta5,exons), tuple(meta6,fasta))
    ch_stats_input = BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_INDEX.out.csi)
        .map { meta, vcf, csi -> [meta, vcf, csi] }
    BCFTOOLS_STATS(ch_stats_input, [[],[]], [[],[]], [[],[]], [[],[]], [[],[]])

    emit:
    vcf      = BCFTOOLS_VIEW.out.vcf       // tuple(meta, vcf.gz)
    index    = BCFTOOLS_INDEX.out.csi       // tuple(meta, vcf.gz.csi)
    stats    = BCFTOOLS_STATS.out.stats     // tuple(meta, stats.txt)
}
