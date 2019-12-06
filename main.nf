#!/usr/bin/env nextflow

params.help= false
params.input_files = false
params.reference = "/projects/data/gatk_bundle/hg19/ucsc.hg19.fasta"						// TODO: remove this hard coded bit
params.output = false
params.skip_split_mnps = false
params.filter = false

def helpMessage() {
    log.info"""
Usage:
    variant_normalization.nf --input_files input_files --reference reference.fasta

This workflow implements a VT VCF normalization pipeline (vt v0.5772)

Input:
    * input_files: the path to a tab-separated values file containing in each row the sample name  and path to the VCF file
    The input file does not have header!
    Example input file:
    sample1	/path/to/your/file.vcf
    sample2	/path/to/your/file2.vcf

Optional input:
    * reference: path to the FASTA genome reference (indexes expected *.fai, *.dict) [default: hg19]
    * output: the folder where to publish output
    * skip_split_mnps: flag indicating not to split MNPs
    * filter: specify the filter to apply if any (e.g.: PASS), only variants with this value will be kept

Output:
    * Normalized VCF file
    * One normalized VCF file per variant type (SNPs, MNPs, indels, BND, other)
    * Tab-separated values file with the absolute paths to the preprocessed BAMs, preprocessed_bams.txt
    * Summary stats and plots on the VCF
    """
}

if (params.help) {
    helpMessage()
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
    cpus 1
    memory '4g'
    //container 'biocontainers/bcftools'
    module 'anaconda/3/2019'
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
    cpus 1
    memory '4g'
    //container 'biocontainers/bcftools'
    module 'anaconda/3/2019'
    tag "${name}"
    publishDir "${publish_dir}/${name}", mode: "copy"

    input:
    	set name, file(vcf) from filtered_vcf

    output:
      set name, file("${vcf.baseName}.normalized.vcf") into normalized_vcf
      set name, file("${vcf.baseName}.normalized.vcf") into normalized_vcf2
      set name, file("${vcf.baseName}.normalized.snps.vcf") into normalized_snps_vcf
      set name, file("${vcf.baseName}.normalized.indels.vcf") into normalized_indels_vcf
      set name, file("${vcf.baseName}.normalized.mnps.vcf") into normalized_mnps_vcf
      set name, file("${vcf.baseName}.normalized.ref.vcf") into normalized_ref_vcf
      set name, file("${vcf.baseName}.normalized.bnd.vcf") into normalized_bnd_vcf
      set name, file("${vcf.baseName}.normalized.other.vcf") into normalized_other_vcf
      file("${vcf.baseName}.normalization.log") into normalization_log

    script:

    """
    # decompose biallelic block substitutions (AC>TG to A>T and C>G)
    # -a: best guess for non blocked substitutions
    # -p: output phased genotypes and PS annotation
    if ["${params.skip_split_mnps}" = true] ; then
      cp $vcf ${vcf.baseName}.atomic.vcf
      touch ${vcf.baseName}.decompose_blocksub_stats.log
    else
      vt decompose_blocksub ${vcf} -a -p -o ${vcf.baseName}.atomic.vcf 2> ${vcf.baseName}.normalization.log
    fi

    # decompose multiallelic variants into biallelic (C>T,G to C>T and C>G)
    vt decompose ${vcf.baseName}.atomic.vcf -o ${vcf.baseName}.biallelic.vcf 2>> ${vcf.baseName}.normalization.log

    # sort the input VCF
    vt sort ${vcf.baseName}.biallelic.vcf -o ${vcf.baseName}.sorted.vcf

    # normalize variants (trim and left alignment)
    vt normalize ${vcf.baseName}.sorted.vcf -r ${params.reference} -o ${vcf.baseName}.trimmed_left_aligned.vcf 2>> ${vcf.baseName}.normalization.log

    # removes duplicated variants
    vt uniq ${vcf.baseName}.trimmed_left_aligned.vcf -o ${vcf.baseName}.normalized.vcf 2>> ${vcf.baseName}.normalization.log

    # separate by variant type once normalized
    bcftools view --types snps -o ${vcf.baseName}.normalized.snps.vcf ${vcf.baseName}.normalized.vcf
    bcftools view --types indels -o ${vcf.baseName}.normalized.indels.vcf ${vcf.baseName}.normalized.vcf
    bcftools view --types mnps -o ${vcf.baseName}.normalized.mnps.vcf ${vcf.baseName}.normalized.vcf
    bcftools view --types ref -o ${vcf.baseName}.normalized.ref.vcf ${vcf.baseName}.normalized.vcf
    bcftools view --types bnd -o ${vcf.baseName}.normalized.bnd.vcf ${vcf.baseName}.normalized.vcf
    bcftools view --types other -o ${vcf.baseName}.normalized.other.vcf ${vcf.baseName}.normalized.vcf

    # delete intermediate files
    rm -f ${vcf.baseName}.atomic.vcf
    rm -f ${vcf.baseName}.biallelic.vcf
    rm -f ${vcf.baseName}.sorted.vcf
    rm -f ${vcf.baseName}.trimmed_left_aligned.vcf
    """
}

process summaryVcf {
  cpus 1
  memory '4g'
  //container 'biocontainers/bcftools'
  module 'anaconda/3/2019'
  tag "${name}"
  publishDir "${publish_dir}/${name}", mode: "copy"

  input:
    set name, file(vcf) from normalized_vcf2

  output:
    file("${vcf.baseName}.stats") into vcf_stats
    file("${name}_stats/*") into vcf_stats_plots

  """
  bcftools stats $vcf > ${vcf.baseName}.stats
  mkdir -p ${name}_stats
  plot-vcfstats -p ${name}_stats --no-PDF --title ${name} ${vcf.baseName}.stats
  """
}

normalized_vcf
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/normalized_vcfs.txt", newLine: true)

normalized_snps_vcf
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/normalized_snps_vcfs.txt", newLine: true)

normalized_indels_vcf
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/normalized_indels_vcfs.txt", newLine: true)

normalized_mnps_vcf
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/normalized_mnps_vcfs.txt", newLine: true)

normalized_bnd_vcf
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/normalized_bnd_vcfs.txt", newLine: true)

normalized_other_vcf
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/normalized_other_vcfs.txt", newLine: true)
