#!/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2

// Define the process
process gatk_index {
        cpus "${params.cpus}"
        debug = true

        // container
        container "${params.gatk__container}"

        input:
        path(params.ref)

        output:
        path("*"), optional:true

        script:
        """
        gatk CreateSequenceDictionary -R ${params.ref} -O ${params.ref}.dict
        """
}
