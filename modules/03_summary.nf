params.cpus = 1
params.memory = "4g"
params.output = "output"


process SUMMARY_VCF {
  cpus params.cpus
  memory params.memory
  tag "${name}"
  publishDir "${params.output}/${name}/metrics", mode: "copy"

  conda (params.enable_conda ? "conda-forge::libgcc-ng=10.3.0 bioconda::bcftools=1.15.1" : null)

  input:
    tuple val(name), file(vcf)

  output:
    tuple val(name), file("${vcf.baseName}.stats*"), emit: vcf_summaries

  """
  bcftools stats $vcf > ${vcf.baseName}.stats
  """
}