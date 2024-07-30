process markDuplicates {
    input:
    tuple val(alpID), path(bam), path(bai)

    output:
    tuple val(alpID), path("${alpID}_q20_rg_d.bam"), path("${alpID}_q20_rg_d.bai"), path("${alpID}_q20_rg_d_metrics.txt")

    script:
    """
    gatk MarkDuplicates \
        CREATE_INDEX=TRUE \
        ASSUME_SORT_ORDER=coordinate \
        QUIET=TRUE \
        I=${bam} \
        O=${alpID}_q20_rg_d.bam \
        M=${alpID}_q20_rg_d_metrics.txt \
        TMP_DIR=temp \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
    samtools index ${alpID}_q20_rg_d.bam ${alpID}_q20_rg_d.bai
    """
}
