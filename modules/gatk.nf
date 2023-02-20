#!/bin/env nextflow

process gatk_dict {
    container 'broadinstitute/gatk:4.3.0.0'

    input:
    path(params.ref)

    output:
    path("*")

    script:
    """
	gatk CreateSequenceDictionary \
		-R ${params.ref} \
		-O ${params.ref}.dict
    """
}
