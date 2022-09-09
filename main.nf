#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2

// Import subworkflows to be run 
include { bwa_index } from './modules/bwa.nf' 
include { samtools_index } from './modules/samtools.nf'
include { gatk_index } from './modules/gatk.nf'

/// Print a header for your pipeline 

log.info """\
     ==============================================
     ==============================================
      I N D E X  R E F E R E N C E  F A S T A - N F  
     ==============================================
     ==============================================
  -._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._      
     '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '. .  
   '.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.    
   : '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '  
   '   '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.   
          `-..,..-'       `-..,..-'       `-..,..-'           
 
                ~~~~ Version: 1.0 ~~~~
 
 Created by the Sydney Informatics Hub, University of Sydney
 Find documentation and more info @ https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf.git
 Cite this pipeline @ INSERT DOI
 Log issues @ https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf/issues
 All default parameters are set in `nextflow.config`.
"""

/// Help function

def helpMessage() {
    log.info"""
  OOPS, YOU FORGOT SOME REQUIRED ARGUMENTS!

  Usage:  nextflow run https://github.com/Sydney-Informatics-Hub/bamQC-nf --ref <reference.fasta> --samtools --gatk --bwa

  Required Arguments:
	
	--ref			Specify path to reference FASTA file 
				(include file extension).
	
  Optional Arguments:
	
	--samtools		Run samtools faidx
	
	--bwa			Run bwa index

	--gatk			Run gatk CreateSequenceDictionary

    """.stripIndent()
}

// run workflow
workflow{

// Show help message if --help is run or if any required params are not 
// provided at runtime

        if ( params.help || params.ref == false ){
        // Invoke the help function above and exit
              helpMessage()
              exit 1

	} else {
	
	if (params.bwa) {
	// create bwa indexes
	bwa_index(params.ref)}

	if (params.samtools) {
	// create samtools index 
	samtools_index(params.ref)}

	if (params.gatk) {
	// create gatk dict file
	gatk_index(params.ref)}
}}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Reference indexes created, runtime info is in `./run_Info`" 
	: "Oops .. something went wrong" )
}
