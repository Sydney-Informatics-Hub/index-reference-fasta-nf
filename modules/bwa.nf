#!/bin/env nextflow

process bwa_index {
	container 'quay.io/biocontainers/bwa:0.7.17--h7132678_9'

	input:
	path(params.ref)

	output:
	path("*")

	script:
	"""
	bwa index \
		-a bwtsw \
		"${params.ref}"
	"""
}
