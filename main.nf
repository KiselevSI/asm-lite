#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Validate required params
if (!params.input) {
    error "ERROR: --input samplesheet.csv is required"
}

if (params.workflow == 'bacteria') {
    include { BACTERIA } from './workflows/bacteria'
} else if (params.workflow == 'phages') {
    include { PHAGES } from './workflows/phages'
} else {
    error "ERROR: --workflow must be 'bacteria' or 'phages', got: '${params.workflow}'"
}

workflow {
    if (params.workflow == 'bacteria') {
        BACTERIA()
    } else {
        PHAGES()
    }
}
