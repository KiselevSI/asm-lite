#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Validate required params
if (!params.input) {
    error "ERROR: --input samplesheet.csv is required"
}
if (!(params.workflow in ['bacteria', 'phages'])) {
    error "ERROR: --workflow must be 'bacteria' or 'phages', got: '${params.workflow}'"
}
if (params.workflow == 'phages' && !params.pharokka_db) {
    error "ERROR: --pharokka_db is required for phages workflow"
}

include { ASMLITE } from './workflows/asmlite'

workflow {
    ASMLITE()
}
