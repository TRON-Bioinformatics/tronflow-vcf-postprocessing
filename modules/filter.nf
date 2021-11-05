params.cpus = 1
params.memory = "4g"
params.output = "."


process FILTER_VCF {
    cpus params.cpus
    memory params.memory
    tag "${name}"


    input:
    	tuple val(name), file(vcf)

    output:
      tuple val(name), file("${vcf.baseName}.filtered.vcf"), emit: filtered_vcfs

    """
    # filter variants
    bcftools view --apply-filters ${params.filter} -o ${vcf.baseName}.filtered.vcf ${vcf}
    """
  }