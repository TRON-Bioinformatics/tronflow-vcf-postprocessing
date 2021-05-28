#!/usr/bin/env nextflow

params.help= false
params.input_files = false
params.reference = false
params.output = false
params.skip_split_mnps = false
params.filter = false
params.decompose_non_blocked_substitutions = false
params.cpus = 1
params.memory = "4g"


if (params.help) {
    log.info params.help
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
      set name, file("${vcf.baseName}.filtered.vcf") into filtered_vcf

    """
    # filter variants
    bcftools view --apply-filter ${params.filter} -o ${vcf.baseName}.filtered.vcf ${vcf}
    """
  }
}
else {
  filtered_vcf = input_files
}

/*
This step sets MAPQ to 0 for all unmapped reads + avoids soft clipping beyond the end of the reference genome
This step reorders chromosomes in the BAM file according to the provided reference (this step is required for GATK)
Adds the required read groups fields to the BAM file. The provided type is added to the BAM sample name.
*/
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
        decomposeNonBlockedSubstitutionsOption = params.decompose_non_blocked_substitutions ? " -a " : ""
        normalizedVcf =  name + ".normalized.vcf"
        decompose_blocksub = params.skip_split_mnps ?
            "" :
            "vt decompose_blocksub - ${decomposeNonBlockedSubstitutionsOption} |"
    """
    # initial sort of the VCF
    vt sort ${vcf} | \
    # decompose biallelic block substitutions (AC>TG to A>T and C>G)
    # -a: best guess for non blocked substitutions
    # -p: output phased genotypes and PS annotation
    ${decompose_blocksub} \
    # decompose multiallelic variants into biallelic (C>T,G to C>T and C>G)
    vt decompose - | \
    # sort the input VCF
    vt sort - | \
    # normalize variants (trim and left alignment)
    vt normalize - -r ${params.reference} | \
    # removes duplicated variants
    vt uniq - -o ${normalizedVcf}
    """
}

process summaryVcf {
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
