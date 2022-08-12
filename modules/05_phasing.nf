params.phasing_sample = false
params.cpus = 1
params.memory = '3g'


process WHATSHAP {
    cpus params.cpus
    memory params.memory
    tag "${patient_name}"
    publishDir "${params.output}/${patient_name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::whatshap=1.4" : null)

    input:
    tuple val(patient_name), file(vcf), val(bams)

    output:
    tuple val(patient_name), file("${patient_name}.phased.vcf"), emit: phased_vcf

    script:
    phasing_bams = bams.toList().stream()
        .filter{ b -> b.split(':')[0] == "${params.phasing_sample}" }
        .collect{ b -> b.split(":")[1] }
        .join(" ") ;
    """
    whatshap phase \
    -o ${patient_name}.phased.vcf \
    --reference=${params.reference} \
    ${vcf} \
    ${phasing_bams}
    """
}