process mapReads {
    // Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = true //turn to false to stop printing command stdout to screen
	tag "WORKING ON: ${params.input}" 
	publishDir "${params.outdir}/processOne", mode: 'symlink'
    container '' 

    input:
    tuple val(alpID), path(fastq1), path(fastq2)

    output:
    tuple val(alpID), path("${alpID}_q20.bam"), path("${alpID}_q20.bai")

    script:
    """
    bowtie2 --local -x ${params.bowtie2_index} \
        -p 12 \
        -1 ${fastq1} \
        -2 ${fastq2} | \
        samtools view -bS - | \
        samtools view -bT ${params.reference} - | \
        samtools sort - -o ${alpID}_q20.bam
    samtools index ${alpID}_q20.bam ${alpID}_q20.bai
    """
}