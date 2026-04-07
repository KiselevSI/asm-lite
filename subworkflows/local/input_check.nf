// subworkflows/local/input_check.nf

workflow INPUT_CHECK {
    take:
    samplesheet // path to CSV

    main:
    reads = channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def fastq_1 = row.fastq_1?.trim()
            def fastq_2 = row.fastq_2?.trim()
            def reference = row.reference?.trim()

            if (!fastq_1) {
                error "ERROR: sample '${row.sample}' is missing fastq_1 in --input"
            }

            if (!params.skip_variants && !reference) {
                error "ERROR: sample '${row.sample}' is missing the 'reference' value in --input. Add a per-sample FASTA path or use --skip_variants."
            }

            def meta = [id: row.sample, single_end: !fastq_2]
            if (reference) {
                def ref_fasta = file(reference, checkIfExists: true)
                meta.reference = ref_fasta.toString()
                meta.ref_name = ref_fasta.baseName
            }

            def reads = fastq_2
                ? [file(fastq_1, checkIfExists: true), file(fastq_2, checkIfExists: true)]
                : [file(fastq_1, checkIfExists: true)]
            [meta, reads]
        }

    emit:
    reads
}
