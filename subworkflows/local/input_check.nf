// subworkflows/local/input_check.nf

workflow INPUT_CHECK {
    take:
    samplesheet // path to CSV

    main:
    reads = channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [id: row.sample, single_end: !row.fastq_2]
            def reads = row.fastq_2
                ? [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
                : [file(row.fastq_1, checkIfExists: true)]
            [meta, reads]
        }

    emit:
    reads
}
