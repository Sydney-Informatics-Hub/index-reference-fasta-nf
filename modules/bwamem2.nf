#!/bin/env nextflow

process bwamem2_index {
	container 'quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0'

	input:
	path(params.ref)

	output:
	path("*")

	script:
	"""
	bwa-mem2 index \
		"${params.ref}"
	"""
}
