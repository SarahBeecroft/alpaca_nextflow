process replaceReadGroups {
    input:
    tuple val(alpID), path(bam), path(bai)

    output:
    tuple val(alpID), path("${alpID}_q20_rg.bam"), path("${alpID}_q20_rg.bai")

    script:
    """
    gatk AddOrReplaceReadGroups \
        QUIET=TRUE \
        I=${bam} \
        O=${alpID}_q20_rg.bam \
        RGID=${alpID} \
        RGLB=${alpID} \
        RGPL=illumina \
        RGSM=${alpID} \
        RGPU=${alpID}
    samtools index ${alpID}_q20_rg.bam ${alpID}_q20_rg.bai
    """
}