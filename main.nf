#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// Print a header for your pipeline 
log.info """\

=======================================================================================
AlpacasAssemble-nf 
=======================================================================================

Created by Sarah Beecroft, Deva Deeptimahanti
#TODO Find documentation @ https://sydney-informatics-hub.github.io/Nextflow_DSL2_template_guide/
#TODO Cite this pipeline @ INSERT DOI

=======================================================================================
Workflow run parameters 
=======================================================================================
input       : ${params.input}
results     : ${params.outDir}
workDir     : ${workflow.workDir}
=======================================================================================

"""

/// Help function 
// This is an example of how to set out the help function that 
// will be run if run command is incorrect or missing. 

def helpMessage() {
    log.info"""
  Usage:  nextflow run main.nf --samples samples.csv --reference path/to/reference.fasta --bowtie2_index_prefix reference_prefix

  Required Arguments:

  --samples		Specify full path and name of sample input file.
  --reference 	Specify full path and name of reference genome fasta file. Assumes fasta.fai index is in the same directory.
  --bowtie2_index	The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.
  --

  Optional Arguments:

  --outDir	Specify path to output directory.
  --bowtie2_indices_path BOWTIE2_INDEXES environment variable--> do I need to include this for nextflow to know where to find the indices

	
""".stripIndent()
}


// Define workflow structure. Include some input/runtime tests here.
// See https://www.nextflow.io/docs/latest/dsl2.html?highlight=workflow#workflow
workflow {

// Show help message if --help is run or (||) a required parameter (input) is not provided

if ( params.help || params.input == false ){   
// Invoke the help function above and exit
	helpMessage()
	exit 1
	// consider adding some extra contigencies here.
	// could validate path of all input files in list?
	// could validate indexes for reference exist?

// If none of the above are a problem, then run the workflow
} else {
	
	// DEFINE CHANNELS 
	// See https://www.nextflow.io/docs/latest/channel.html#channels
	// See https://training.nextflow.io/basic_training/channels/ 
	Channel
    .fromPath(params.samples)
    .splitCsv(header:true)
    .map { row -> tuple(row.alpID, file(row.fastq1), file(row.fastq2)) }
    .set { samples_ch }

	// VALIDATE INPUT SAMPLES 
	check_input(Channel.fromPath(params.input, checkIfExists: true))

	// QC AND TRIM FASTQ
	qcAndTrim(samples_ch)
	
	// MAP READS USING BOWTIE2
    mapReads(qcAndTrim.out)
	
	// REPLACE READGROUP INFO
    replaceReadGroups(mapReads.out)
	
	// MARK DUPLICATES WITH GATK
    markDuplicates(replaceReadGroups.out)
	
	// CALL VARIANTS WITH GATK
    callVariants(markDuplicates.out)

}}

// Print workflow execution summary 
workflow.onComplete {
summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
results     : ${params.outDir}

=======================================================================================
  """
println summary

}

