params.cpus = 1
params.memory = "4g"
params.output = "output"
params.filter = false


process FIX_VCF_HEADER {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.0" : null)

    input:
      tuple val(name), file(vcf)

    output:
      tuple val(name), file("${name}.reheader.vcf"), emit: reheader_vcfs

    script:
    header = params.header ? "-H ${params.header}" : ""
    """
    gatk FixVcfHeader ${header} -I ${vcf} -O ${name}.reheader.vcf
    """
}
