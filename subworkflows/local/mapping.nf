// subworkflows/local/mapping.nf

include { BWA_INDEX as BWA_INDEX_REF           } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_MEM_REF               } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_REF  } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX                        } from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_MARKDUPLICATES                 } from '../../modules/nf-core/picard/markduplicates/main'

workflow MAPPING {
    take:
    reads       // tuple(meta[+reference,+ref_name], [fastq_1, fastq_2])

    main:

    ch_reads_by_ref = reads
        .map { meta, reads_files ->
            if (!meta.reference) {
                error "ERROR: sample '${meta.id}' does not have a reference in --input. Add column 'reference' or use --skip_variants."
            }
            [meta.reference, meta, reads_files]
        }

    // Unique reference FASTA files declared in the samplesheet
    ch_refs = reads
        .map { meta, reads_files ->
            if (!meta.reference) {
                error "ERROR: sample '${meta.id}' does not have a reference in --input. Add column 'reference' or use --skip_variants."
            }
            [meta.reference, meta.ref_name]
        }
        .unique()
        .map { ref_path, ref_name ->
            def ref_meta = [id: ref_name, ref_key: ref_path]
            [ref_meta, file(ref_path, checkIfExists: true)]
        }

    // Index each reference
    BWA_INDEX_REF(ch_refs)

    // Create faidx for each reference (needed by picard)
    ch_faidx_input = ch_refs.map { ref_meta, fasta -> [ref_meta, fasta, []] }
    SAMTOOLS_FAIDX(ch_faidx_input, false)

    ch_ref_resources = ch_refs
        .map { ref_meta, fasta -> [ref_meta.ref_key, ref_meta, fasta] }
        .combine(
            BWA_INDEX_REF.out.index.map { ref_meta, index -> [ref_meta.ref_key, index] },
            by: 0
        )
        .combine(
            SAMTOOLS_FAIDX.out.fai.map { ref_meta, fai -> [ref_meta.ref_key, fai] },
            by: 0
        )
        .map { ref_key, ref_meta, fasta, index, fai ->
            [ref_key, ref_meta, fasta, index, fai]
        }

    ch_for_bwa = ch_reads_by_ref
        .combine(ch_ref_resources, by: 0)
        .map { ref_key, meta, reads_files, ref_meta, ref_fasta, index, fai ->
            [meta, reads_files, ref_meta, index, ref_fasta]
        }

    BWA_MEM_REF(
        ch_for_bwa.map { meta, reads_files, ref_meta, index, ref_fasta -> [meta, reads_files] },
        ch_for_bwa.map { meta, reads_files, ref_meta, index, ref_fasta -> [ref_meta, index] },
        ch_for_bwa.map { meta, reads_files, ref_meta, index, ref_fasta -> [ref_meta, ref_fasta] },
        true  // sort BAM
    )

    // Mark duplicates
    ch_ref_fasta_fai = ch_ref_resources
        .map { ref_key, ref_meta, fasta, index, fai -> [ref_key, ref_meta, fasta, fai] }

    ch_dedup_with_ref = BWA_MEM_REF.out.bam
        .map { meta, bam -> [meta.reference, meta, bam] }
        .combine(ch_ref_fasta_fai, by: 0)
        .map { ref_key, meta, bam, ref_meta, fasta, fai ->
            [meta, bam, ref_meta, fasta, fai]
        }

    PICARD_MARKDUPLICATES(
        ch_dedup_with_ref.map { meta, bam, ref_meta, fasta, fai -> [meta, bam] },
        ch_dedup_with_ref.map { meta, bam, ref_meta, fasta, fai -> [ref_meta, fasta, fai] }
    )

    // Index deduped BAM
    SAMTOOLS_INDEX_REF(PICARD_MARKDUPLICATES.out.bam)

    // Join BAM + BAI
    ch_bam_bai = PICARD_MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX_REF.out.index)
        .map { meta, bam, bai -> [meta, bam, bai] }

    // BAM + BAI + ref fasta + fai for variant calling
    ch_bam_bai_ref = ch_bam_bai
        .map { meta, bam, bai -> [meta.reference, meta, bam, bai] }
        .combine(ch_ref_fasta_fai, by: 0)
        .map { ref_key, meta, bam, bai, ref_meta, fasta, fai -> [meta, bam, bai, fasta, fai] }

    emit:
    bam_bai     = ch_bam_bai        // tuple(meta[+ref_name], bam, bai)
    bam_bai_ref = ch_bam_bai_ref    // tuple(meta[+ref_name], bam, bai, ref_fasta, ref_fai)
}
