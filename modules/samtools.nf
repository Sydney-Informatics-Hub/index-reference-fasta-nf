#!/bin/env nextflow

process samtools_index {
	container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
	
	input:
	path(params.ref)

	output:
	path("*")

	script:
	"""
	samtools faidx \
		"${params.ref}"
	"""
}
