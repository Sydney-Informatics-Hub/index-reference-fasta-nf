#!/bin/env nextflow

process gatk_dict {
        container "${params.gatk__container}"

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
