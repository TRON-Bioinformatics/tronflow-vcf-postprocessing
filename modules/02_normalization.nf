params.cpus = 1
params.memory = "4g"
params.output = "output"
params.vcf_without_ad = false


process BCFTOOLS_NORM {
    cpus params.cpus
    memory params.memory
    tag "${name}"

    conda (params.enable_conda ? "conda-forge::libgcc-ng=10.3.0 conda-forge::gsl=2.7 bioconda::bcftools=1.15.1" : null)

    input:
    	tuple val(name), file(vcf)

    output:
      tuple val(name), file("${name}.atomic_left_aligned.vcf"), emit: normalized_vcfs

    script:
    keep_ad_sum = params.vcf_without_ad ? "--keep-sum AD" : ""
    """
    # initial sort of the VCF
    bcftools sort ${vcf} | \

    # checks reference genome, decompose multiallelics, trim and left align indels
    bcftools norm --multiallelics -any ${keep_ad_sum} --check-ref e --fasta-ref ${params.reference} \
    --old-rec-tag OLD_CLUMPED - > ${name}.atomic_left_aligned.vcf
    """
}

process VT_DECOMPOSE_COMPLEX {
    cpus params.cpus
    memory params.memory
    tag "${name}"

    conda (params.enable_conda ? "openssl=3.0.0 bioconda::vt=0.57721" : null)

    input:
    	tuple val(name), file(vcf)

    output:
      tuple val(name), file("${name}.decomposed.vcf"), emit: decomposed_vcfs

    script:
    """
    vt decompose_blocksub -a -p ${vcf} > ${name}.decomposed.vcf
    """
}

process REMOVE_DUPLICATES {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda (params.enable_conda ? "conda-forge::libgcc-ng=10.3.0 conda-forge::gsl=2.7 bioconda::bcftools=1.15.1" : null)

    input:
    	tuple val(name), file(vcf)

    output:
      tuple val(name), file("${name}.normalized.vcf"), emit: deduplicated_vcfs

    script:
    """
    # remove duplicates after normalisation
    bcftools norm --rm-dup exact -o ${name}.normalized.vcf ${vcf}
    """
}