process callVariants {
    input:
    tuple val(alpID), path(bam), path(bai), path(metrics)

    output:
    path("${alpID}.g.vcf")

    script:
    """
    gatk HaplotypeCaller \
        -R ${params.reference} \
        -I ${bam} \
        -ERC GVCF \
        -O ${alpID}.g.vcf
    """
}