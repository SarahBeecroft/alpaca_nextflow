process qcAndTrim {

	debug = true //turn to false to stop printing command stdout to screen
	tag "WORKING ON: ${params.input}" 
	publishDir "${params.outdir}/processOne", mode: 'symlink'
  	container '' 


    input:
    tuple val(alpID), path(fastq1), path(fastq2)

    output:
    tuple val(alpID), path("${alpID}_1.fastq"), path("${alpID}_2.fastq"), path("${alpID}_fastp.html")

    script:
    """
    fastp --thread 8 \
        -h report/${alpID}_fastp.html \
        -i ${fastq1} \
        -I ${fastq2} \
        -o qc_reads/${alpID}_1.fastq \
        -O qc_reads/${alpID}_2.fastq
    """
}