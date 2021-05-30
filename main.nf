#!/usr/bin/env nextflow

params.help= false
params.input_files = false
params.reference = false
params.output = false
params.skip_decompose_complex = false
params.filter = false
params.cpus = 1
params.memory = "4g"


if (params.help) {
    log.info params.help_message
    exit 0
}

publish_dir = "output"
if (params.output) {
  publish_dir = params.output
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'vcf'], sep: "\t")
    .map{ row-> tuple(row.name, file(row.vcf)) }
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}

if (params.filter) {
  process filterVcf {
    cpus params.cpus
    memory params.memory
    tag "${name}"


    input:
    	set name, file(vcf) from input_files

    output:
      set name, file("${vcf.baseName}.filtered.vcf") into filtered_vcf, filtered_vcf_for_stats

    """
    # filter variants
    bcftools view --apply-filter ${params.filter} -o ${vcf.baseName}.filtered.vcf ${vcf}
    """
  }
}
else {
    input_files.into { filtered_vcf; filtered_vcf_for_stats }
}

process summaryVcfBefore {
  cpus params.cpus
  memory params.memory
  tag "${name}"
  publishDir "${publish_dir}/${name}/metrics", mode: "copy"

  input:
    set name, file(vcf) from filtered_vcf_for_stats

  output:
    file("${vcf.baseName}.stats*")

  """
  bcftools stats $vcf > ${vcf.baseName}.stats
  """
}

process normalizeVcf {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${publish_dir}/${name}", mode: "copy"

    input:
    	set name, file(vcf) from filtered_vcf

    output:
      set name, file("${name}.normalized.vcf") into normalized_vcf2
      set name, val("${publish_dir}/${name}/${name}.normalized.vcf") into normalized_vcf

    script:
        //decompose_complex = params.skip_decompose_complex ? "" : "bcftools norm --atomize - |"
        decompose_complex = params.skip_decompose_complex ? "" : "vt decompose_blocksub -a -p - |"

    """
    # initial sort of the VCF
    bcftools sort ${vcf} | \

    # checks reference genome, decompose multiallelics, trim and left align indels
    bcftools norm --multiallelics -any --keep-sum AD --check-ref e --fasta-ref ${params.reference} \
    --old-rec-tag OLD_CLUMPED - | \

    # decompose complex variants
    ${decompose_complex}

    # remove duplicates after normalisation
    bcftools norm --rm-dup exact -o ${name}.normalized.vcf -
    """
}

process summaryVcfAfter {
  cpus params.cpus
  memory params.memory
  tag "${name}"
  publishDir "${publish_dir}/${name}/metrics", mode: "copy"

  input:
    set name, file(vcf) from normalized_vcf2

  output:
    file("${vcf.baseName}.stats*")

  """
  bcftools stats $vcf > ${vcf.baseName}.stats
  """
}

normalized_vcf
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/normalized_vcfs.txt", newLine: true)
