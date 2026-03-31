// subworkflows/local/mapping.nf

include { BWA_INDEX as BWA_INDEX_REF           } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_MEM_REF               } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_REF  } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX                        } from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_MARKDUPLICATES                 } from '../../modules/nf-core/picard/markduplicates/main'

workflow MAPPING {
    take:
    reads       // tuple(meta, [fastq_1, fastq_2])
    references  // path to references.csv

    main:

    // Parse references CSV
    ch_refs = channel
        .fromPath(references)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def ref_meta = [id: row.name]
            [ref_meta, file(row.fasta, checkIfExists: true)]
        }

    // Index each reference
    BWA_INDEX_REF(ch_refs)

    // Create faidx for each reference (needed by picard)
    ch_faidx_input = ch_refs.map { ref_meta, fasta -> [ref_meta, fasta, []] }
    SAMTOOLS_FAIDX(ch_faidx_input, false)

    // Combine: each sample x each reference
    ch_reads_x_refs = reads
        .combine(ch_refs)
        .map { meta, reads_files, ref_meta, ref_fasta ->
            def new_meta = meta + [ref_name: ref_meta.id]
            [new_meta, reads_files, ref_meta, ref_fasta]
        }

    // Join with BWA index by ref_meta.id
    ch_for_bwa = ch_reads_x_refs
        .map { meta, reads_files, ref_meta, ref_fasta ->
            [ref_meta.id, meta, reads_files, ref_meta, ref_fasta]
        }
        .combine(
            BWA_INDEX_REF.out.index.map { ref_meta, index -> [ref_meta.id, index] },
            by: 0
        )
        .map { ref_id, meta, reads_files, ref_meta, ref_fasta, index ->
            [meta, reads_files, ref_meta, index, ref_fasta]
        }

    BWA_MEM_REF(
        ch_for_bwa.map { meta, reads_files, ref_meta, index, ref_fasta -> [meta, reads_files] },
        ch_for_bwa.map { meta, reads_files, ref_meta, index, ref_fasta -> [ref_meta, index] },
        ch_for_bwa.map { meta, reads_files, ref_meta, index, ref_fasta -> [ref_meta, ref_fasta] },
        true  // sort BAM
    )

    // Mark duplicates
    // Get fai from SAMTOOLS_FAIDX
    ch_ref_fasta_fai = ch_refs
        .join(SAMTOOLS_FAIDX.out.fai)
        .map { ref_meta, fasta, fai -> [ref_meta.id, ref_meta, fasta, fai] }

    ch_dedup_with_ref = BWA_MEM_REF.out.bam
        .map { meta, bam -> [meta.ref_name, meta, bam] }
        .combine(ch_ref_fasta_fai, by: 0)
        .map { ref_id, meta, bam, ref_meta, fasta, fai ->
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
        .map { meta, bam, bai -> [meta.ref_name, meta, bam, bai] }
        .combine(ch_refs.map { ref_meta, fasta -> [ref_meta.id, fasta] }, by: 0)
        .combine(SAMTOOLS_FAIDX.out.fai.map { ref_meta, fai -> [ref_meta.id, fai] }, by: 0)
        .map { ref_id, meta, bam, bai, fasta, fai -> [meta, bam, bai, fasta, fai] }

    emit:
    bam_bai     = ch_bam_bai        // tuple(meta[+ref_name], bam, bai)
    bam_bai_ref = ch_bam_bai_ref    // tuple(meta[+ref_name], bam, bai, ref_fasta, ref_fai)
}
