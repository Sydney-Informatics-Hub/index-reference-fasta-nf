#!/bin/env nextflow 

// Enable DSL-2 syntax
nextflow.enable.dsl=2

// Define the process
process bwa_index {
        cpus "${params.cpus}"
	debug = true

	// container
	container "${params.bwa__container}"

	input:
	path(params.ref)	

	output:
	path("*"), optional:true
	//tuple val(meta), path("*.bai") , optional:true, emit: bai

	script:
	"""
	bwa index -a bwtsw "${params.ref}"
	"""
}
