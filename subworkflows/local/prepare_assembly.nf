// subworkflows/local/prepare_assembly.nf

include { UNICYCLER      } from '../../modules/nf-core/unicycler/main'
include { GUNZIP         } from '../../modules/nf-core/gunzip/main'
include { BWA_INDEX      } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM        } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow PREPARE_ASSEMBLY {
    take:
    reads   // tuple(meta, [fastq_1, fastq_2])

    main:
    ch_reads_by_sample = reads
        .map { meta, fastqs -> [meta.id, meta, fastqs] }

    // 1. Assemble with Unicycler
    // Input: tuple(meta, shortreads, longreads) — longreads=[] for none
    ch_unicycler_input = reads.map { meta, fastqs -> [meta, fastqs, []] }
    UNICYCLER(ch_unicycler_input)

    // 2. Gunzip assembly (unicycler outputs .fa.gz)
    GUNZIP(UNICYCLER.out.scaffolds)
    ch_assembly_by_sample = GUNZIP.out.gunzip
        .map { meta, fasta -> [meta.id, meta, fasta] }

    // 3. Index assembly for BWA
    BWA_INDEX(GUNZIP.out.gunzip)
    ch_index_by_sample = BWA_INDEX.out.index
        .map { meta, index -> [meta.id, meta, index] }

    // 4. Map reads back to assembly
    // BWA_MEM(tuple(meta,reads), tuple(meta2,index), tuple(meta3,fasta), sort_bam)
    ch_bwa_input = ch_reads_by_sample
        .join(ch_assembly_by_sample, by: 0)
        .join(ch_index_by_sample, by: 0)
        .map { sample_id, meta, fastqs, asm_meta, fasta, idx_meta, index ->
            [meta, fastqs, asm_meta, fasta, idx_meta, index]
        }

    BWA_MEM(
        ch_bwa_input.map { meta, fastqs, asm_meta, fasta, idx_meta, index -> [meta, fastqs] },
        ch_bwa_input.map { meta, fastqs, asm_meta, fasta, idx_meta, index -> [idx_meta, index] },
        ch_bwa_input.map { meta, fastqs, asm_meta, fasta, idx_meta, index -> [asm_meta, fasta] },
        true    // sort BAM
    )

    // 5. Index BAM
    SAMTOOLS_INDEX(BWA_MEM.out.bam)

    // 6. Join BAM + BAI
    ch_bam_bai = BWA_MEM.out.bam
        .join(SAMTOOLS_INDEX.out.index)
        .map { meta, bam, bai -> [meta, bam, bai] }

    emit:
    scaffolds = GUNZIP.out.gunzip   // tuple(meta, fasta) — unzipped
    bam_bai   = ch_bam_bai          // tuple(meta, bam, bai)
}
