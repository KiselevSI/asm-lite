#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Validate required params
if (!params.input) {
    error "ERROR: --input samplesheet.csv is required"
}
if (params.pharokka && !params.pharokka_db) {
    error "ERROR: --pharokka_db is required when --pharokka is enabled"
}

include { ASMLITE } from './workflows/asmlite'

workflow {
    ASMLITE()
}
