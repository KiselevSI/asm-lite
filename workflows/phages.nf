// workflows/phages.nf

include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { FASTP             } from '../modules/nf-core/fastp/main'
include { FASTQC            } from '../modules/nf-core/fastqc/main'
include { KRAKEN2_KRAKEN2   } from '../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN   } from '../modules/nf-core/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS } from '../modules/nf-core/bracken/combinebrackenoutputs/main'
include { BRACKEN_ADD_UNCLASSIFIED      } from '../modules/local/bracken_add_unclassified'
include { PREPARE_ASSEMBLY  } from '../subworkflows/local/prepare_assembly'
include { PHAROKKA_PHAROKKA } from '../modules/nf-core/pharokka/pharokka/main'
include { MAPPING           } from '../subworkflows/local/mapping'
include { CALL_VARIANTS     } from '../subworkflows/local/call_variants'
include { STATS             } from '../subworkflows/local/stats'
include { MULTIQC           } from '../modules/nf-core/multiqc/main'

workflow PHAGES {
    main:
    ch_multiqc     = channel.empty()

    // 1. Parse samplesheet
    INPUT_CHECK(params.input)
    ch_reads = INPUT_CHECK.out.reads

    // 2. Trim
    ch_fastp_input = ch_reads.map { meta, reads -> [meta, reads, []] }
    FASTP(ch_fastp_input, false, false, false)
    ch_trimmed = FASTP.out.reads
    ch_multiqc  = ch_multiqc.mix(FASTP.out.json.map { meta, json -> json })

    // 3. QC
    if (!params.skip_qc) {
        FASTQC(ch_trimmed)
        ch_multiqc  = ch_multiqc.mix(FASTQC.out.zip.map { meta, zip -> zip })
    }

    // 4. Kraken
    if (!params.skip_kraken && params.kraken2_db) {
        ch_kraken_db = channel.value(file(params.kraken2_db, checkIfExists: true))

        KRAKEN2_KRAKEN2(ch_trimmed, ch_kraken_db, false, false)

        BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report, ch_kraken_db)

        ch_bracken_kraken = BRACKEN_BRACKEN.out.reports
            .join(KRAKEN2_KRAKEN2.out.report)
        BRACKEN_ADD_UNCLASSIFIED(ch_bracken_kraken)

        ch_all_bracken = BRACKEN_BRACKEN.out.reports.map { meta, report -> report }.collect()
        BRACKEN_COMBINEBRACKENOUTPUTS(
            [[id: 'all_samples'], ch_all_bracken]
        )
    }

    // 5. Assembly + map to self
    PREPARE_ASSEMBLY(ch_trimmed)

    // 6. Annotate with Pharokka
    // PHAROKKA_PHAROKKA(tuple(meta, fasta), pharokka_db)
    ch_pharokka_db = channel.value(file(params.pharokka_db, checkIfExists: true))
    PHAROKKA_PHAROKKA(PREPARE_ASSEMBLY.out.scaffolds, ch_pharokka_db)

    // 7. Stats on assembly
    STATS(PREPARE_ASSEMBLY.out.bam_bai, PREPARE_ASSEMBLY.out.scaffolds)
    ch_multiqc  = ch_multiqc.mix(STATS.out.wgs_metrics.map { meta, f -> f })
    ch_multiqc  = ch_multiqc.mix(STATS.out.alignment_metrics.map { meta, f -> f })
    ch_multiqc  = ch_multiqc.mix(STATS.out.samtools_stats.map { meta, f -> f })
    ch_multiqc  = ch_multiqc.mix(STATS.out.samtools_flagstat.map { meta, f -> f })
    ch_multiqc  = ch_multiqc.mix(STATS.out.quast_results.map { meta, d -> d })

    // 8. Map to references + variant calling
    if (!params.skip_variants && params.references) {
        MAPPING(ch_trimmed, params.references)

        CALL_VARIANTS(MAPPING.out.bam_bai_ref)
        ch_multiqc  = ch_multiqc.mix(CALL_VARIANTS.out.stats.map { meta, f -> f })
    }

    // 9. MultiQC
    ch_multiqc_config = params.multiqc_config ? channel.fromPath(params.multiqc_config) : channel.empty()
    MULTIQC(
        ch_multiqc.collect().map { files ->
            [[id: 'multiqc'], files, [], [], [], []]
        }
    )
}
