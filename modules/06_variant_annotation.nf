params.snpeff_memory = "3g"
params.memory = "3g"
params.cpus = 1
params.output = "."
params.snpeff_datadir = false
params.snpeff_organism = false
params.snpeff_args = ""


process VARIANT_ANNOTATION_SNPEFF {
    cpus params.cpus
    memory params.snpeff_memory
    publishDir "${params.output}/${name}", mode: "copy"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::snpeff=5.0" : null)

    input:
        tuple val(name), file(vcf)

    output:
    tuple val(name), file("${name}.annotated.vcf") , emit: annotated_vcf

    script:
    datadir_arg = params.snpeff_datadir ? "-dataDir ${params.snpeff_datadir}" : ""
    """
    snpEff -Xmx${params.snpeff_memory} eff \
    ${datadir_arg} ${params.snpeff_args} \
    -nodownload ${params.snpeff_organism} ${vcf} > ${name}.annotated.vcf
    """
}

process VARIANT_ANNOTATION_BCFTOOLS {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/${name}", mode: "copy"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::bcftools=1.15.1" : null)

    input:
        tuple val(name), file(vcf)

    output:
    tuple val(name), file("${name}.annotated.vcf") , emit: annotated_vcf

    """
    bcftools csq \\
        --fasta-ref ${params.reference} \\
        --gff-annot ${params.gff} ${vcf} \\
        --output-type v \\
        --output ${name}.annotated.vcf \\
        --phase a
    """
}
