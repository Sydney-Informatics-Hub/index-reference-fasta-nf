#!/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2

// Define the process
process samtools_index {
        cpus "${params.cpus}"
        debug = true

        // container
        container "${params.samtools__container}"

        input:
        path(params.ref)

        output:
        path("*"), optional:true

        script:
        """
        samtools faidx "${params.ref}"
        """
}
