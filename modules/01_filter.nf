params.cpus = 1
params.memory = "4g"
params.output = "output"
params.filter = false


process FILTER_VCF {
    cpus params.cpus
    memory params.memory
    tag "${name}"

    conda (params.enable_conda ? "conda-forge::libgcc-ng=10.3.0 conda-forge::gsl=2.7 bioconda::bcftools=1.15.1" : null)

    input:
    	tuple val(name), file(vcf)

    output:
      tuple val(name), file("${vcf.baseName}.filtered.vcf"), emit: filtered_vcfs

    """
    # filter variants
    bcftools view --apply-filters ${params.filter} -o ${vcf.baseName}.filtered.vcf ${vcf}
    """
  }