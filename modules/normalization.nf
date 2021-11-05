params.cpus = 1
params.memory = "4g"
params.output = "output"
params.skip_decompose_complex = false
params.vcf_without_ad = false


process NORMALIZE_VCF {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    input:
    	tuple val(name), file(vcf)

    output:
      tuple val(name), file("${name}.normalized.vcf"), emit: normalized_vcfs

    script:
        //decompose_complex = params.skip_decompose_complex ? "" : "bcftools norm --atomize - |"
        decompose_complex = params.skip_decompose_complex ? "" : "vt decompose_blocksub -a -p - |"
        keep_ad_sum = params.vcf_without_ad ? "--keep-sum AD" : ""
    """
    # initial sort of the VCF
    bcftools sort ${vcf} | \

    # checks reference genome, decompose multiallelics, trim and left align indels
    bcftools norm --multiallelics -any ${keep_ad_sum} --check-ref e --fasta-ref ${params.reference} \
    --old-rec-tag OLD_CLUMPED - | \

    # decompose complex variants
    #${decompose_complex}

    # remove duplicates after normalisation
    bcftools norm --rm-dup exact -o ${name}.normalized.vcf -
    """
}