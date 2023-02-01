params.cpus = 1
params.memory = "4g"
params.output = ""
params.mapping_quality = 0
params.base_call_quality = 0
params.enable_conda = false


process VAFATOR {
    cpus params.cpus
    memory params.memory
    tag "${patient_name}"
    publishDir "${params.output}/${patient_name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::vafator=2.2.0" : null)

    input:
    tuple val(patient_name), file(vcf), val(bams), val(purities), val(clonalities)

    output:
    tuple val(patient_name), file("${patient_name}.vaf.vcf"), emit: annotated_vcf

    script:
    bams_param = bams.collect { b -> "--bam " + b.split(":").join(" ") }.join(" ")
    purity_param = purities.collect { b -> "--purity " + b.split(":").join(" ") }.join(" ")
    clonality_param = clonalities.collect { b -> "--tumor-ploidy " + b.split(":").join(" ") }.join(" ")
    """
    vafator \
    --input-vcf ${vcf} \
    --output-vcf ${patient_name}.vaf.vcf \
    --mapping-quality ${params.mapping_quality} \
    --base-call-quality ${params.base_call_quality} \
    ${bams_param} ${purity_param} ${clonality_param}
    """
}


process MULTIALLELIC_FILTER {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::vafator=2.2.0" : null)

    input:
    tuple val(name), file(vcf)

    output:
    tuple val(name), file("${name}.filtered_multiallelics.vcf"), emit: filtered_vcf

    script:
    """
    multiallelics-filter --input-vcf ${vcf} --output-vcf ${name}.filtered_multiallelics.vcf
    """
}
