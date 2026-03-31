// subworkflows/local/stats.nf

include { PICARD_COLLECTWGSMETRICS              } from '../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../../modules/nf-core/picard/collectalignmentsummarymetrics/main'
include { SAMTOOLS_STATS                        } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT                     } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_ASM  } from '../../modules/nf-core/samtools/faidx/main'
include { QUAST                                 } from '../../modules/nf-core/quast/main'
include { CONTIG_ABUNDANCE                      } from '../../modules/local/contig_abundance'

workflow STATS {
    take:
    bam_bai     // tuple(meta, bam, bai) — reads mapped to own assembly
    assembly    // tuple(meta, fasta)     — assembly

    main:
    // Create faidx for assembly (needed by picard and samtools)
    ch_faidx_input = assembly.map { meta, fasta -> [meta, fasta, []] }
    SAMTOOLS_FAIDX_ASM(ch_faidx_input, false)

    // picard collectwgsmetrics: tuple(meta, bam, bai), tuple(meta2, fasta), tuple(meta3, fai), intervallist
    PICARD_COLLECTWGSMETRICS(
        bam_bai,
        assembly,
        SAMTOOLS_FAIDX_ASM.out.fai,
        []
    )

    // picard collectalignmentsummarymetrics: tuple(meta, bam), tuple(meta2, fasta)
    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(
        bam_bai.map { meta, bam, bai -> [meta, bam] },
        assembly
    )

    // samtools stats: tuple(meta, input, input_index), tuple(meta2, fasta, fai)
    ch_fasta_fai = assembly
        .join(SAMTOOLS_FAIDX_ASM.out.fai)
        .map { meta, fasta, fai -> [meta, fasta, fai] }

    SAMTOOLS_STATS(bam_bai, ch_fasta_fai)

    // samtools flagstat: tuple(meta, bam, bai)
    SAMTOOLS_FLAGSTAT(bam_bai)

    // QUAST: tuple(meta, fasta), tuple(meta2, fasta_ref), tuple(meta3, gff)
    QUAST(assembly, [[],[]], [[],[]])

    // Contig abundance: tuple(meta, bam, bai)
    CONTIG_ABUNDANCE(bam_bai)

    emit:
    wgs_metrics       = PICARD_COLLECTWGSMETRICS.out.metrics              // tuple(meta, metrics)
    alignment_metrics = PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics // tuple(meta, metrics)
    samtools_stats    = SAMTOOLS_STATS.out.stats                          // tuple(meta, stats)
    samtools_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat                    // tuple(meta, flagstat)
    quast_results     = QUAST.out.results                                // tuple(meta, dir)
    contig_abundance  = CONTIG_ABUNDANCE.out.abundance                   // tuple(meta, tsv)
}
