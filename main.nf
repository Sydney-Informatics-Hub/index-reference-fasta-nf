#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// subworkflows to be run
include { bwa_index	} from './modules/bwa.nf'
include { bwamem2_index	} from './modules/bwamem2.nf'
include { samtools_index} from './modules/samtools.nf'
include { gatk_dict	} from './modules/gatk.nf'

// Print header
log.info """\

===================================================================
I N D E X  R E F E R E N C E  F A S T A - N F
===================================================================
Created by the Sydney Informatics Hub, University of Sydney
Documentation   @ https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf
Log issues      @ https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf/issues
===================================================================
Workflow run parameters
===================================================================
reference: ${params.ref}
workDir: ${workflow.workDir}
whoami: ${params.whoami}
gadi_account: ${params.gadi_account}
bwa: ${params.bwa}
bwamem2: ${params.bwamem2}
samtools: ${params.samtools}
gatk: ${params.gatk}
===================================================================
"""

// Help function
def helpMessage() {
    log.info"""
OOPS, YOU FORGOT SOME REQUIRED ARGUMENTS!

Usage:  nextflow run main.nf --ref <reference.fasta> --gatk --bwa --bwamem2 --samtools

Required Arguments:

	--ref       Specify path to reference FASTA file
                (include file extension).

Optional Arguments:

	--bwa       Run bwa index

	--bwamem2 Run bwa-mem2 index

	--gatk      Run gatk CreateSequenceDictionary

	--samtools  Run samtools faidx
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

	if (params.samtools) {
	// create samtools index
	samtools_index(params.ref)
	}

	if (params.bwa) {
	// create bwa indexes
	bwa_index(params.ref)}

	if (params.bwamem2) {
	// create bwa indexes
	bwamem2_index(params.ref)}

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
