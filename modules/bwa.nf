#!/bin/env nextflow 

process bwa_index {
	container "${params.bwa__container}"

	input:
	path(params.ref)	

	output:
	path("*"), optional:true
	//tuple val(meta), path("*.bai") , optional:true, emit: bai

	script:
	"""
	bwa index \
		-a bwtsw \
		"${params.ref}"
	"""
}
