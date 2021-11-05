params.cpus = 1
params.memory = "4g"
params.output = "."


process SUMMARY_VCF {
  cpus params.cpus
  memory params.memory
  tag "${name}"
  publishDir "${params.output}/${name}/metrics", mode: "copy"

  input:
    tuple val(name), file(vcf)

  output:
    tuple val(name), file("${vcf.baseName}.stats*"), emit: vcf_summaries

  """
  bcftools stats $vcf > ${vcf.baseName}.stats
  """
}