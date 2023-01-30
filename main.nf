#!/usr/bin/env nextflow


nextflow.enable.dsl=2

// subworkflows to be run 
include { bwa_index	} from './modules/bwa.nf' 
include { samtools_index} from './modules/samtools.nf'
include { gatk_dict	} from './modules/gatk.nf'

// Print header

log.info """\

===================================================================
I N D E X  R E F E R E N C E  F A S T A - N F                   
===================================================================
Created by the Sydney Informatics Hub, University of Sydney
Documentation		@ https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf
Log issues		@ https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf/issues
===================================================================
Workflow run parameters 
===================================================================
version		: ${params.version}
reference	: ${params.ref}
workDir		: ${workflow.workDir}
===================================================================
"""

// Help function

def helpMessage() {
    log.info"""
  OOPS, YOU FORGOT SOME REQUIRED ARGUMENTS!

  Usage:  nextflow run https://github.com/Sydney-Informatics-Hub/bamQC-nf --ref <reference.fasta> --samtools --gatk --bwa

  Required Arguments:
	
	--ref			Specify path to reference FASTA file 
				(include file extension).
	
  Optional Arguments:
	
	--bwa			Run bwa index

	--gatk			Run gatk CreateSequenceDictionary

    """.stripIndent()
}

// run workflow
workflow{

// Help message

if ( params.help || params.ref == false ){
	// Invoke the help function above and exit
	helpMessage()
	exit 1

} else {

	// create samtools index (default)
	samtools_index(params.ref)
	
	if (params.bwa) {
	// create bwa indexes
	bwa_index(params.ref)}

	if (params.gatk) {
	// create gatk dictionary
	gatk_dict(params.ref)}
}}

workflow.onComplete {

summary = """
===================================================================
Workflow execution summary
===================================================================
Duration	: ${workflow.duration}
Success		: ${workflow.success}
workDir		: ${workflow.workDir}
Exit status	: ${workflow.exitStatus}
===================================================================
"""

println summary

}
