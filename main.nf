#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NORMALIZE_VCF } from './modules/normalization'
include { FILTER_VCF } from './modules/filter'
include { SUMMARY_VCF; SUMMARY_VCF as SUMMARY_VCF_2 } from './modules/summary'

params.help= false
params.input_files = false
params.input_vcf = false
params.reference = false
params.output = "output"
params.skip_decompose_complex = false
params.filter = false
params.cpus = 1
params.memory = "4g"
params.vcf_without_ad = false


if (params.help) {
    log.info params.help_message
    exit 0
}

if (! params.input_files && ! params.input_vcf) {
  exit 1, "Neither --input-files or --input-vcf are provided!"
}
else if (params.input_files && params.input_vcf) {
  exit 1, "Both --input-files and --input-vcf are provided! Please, provide only one."
}
else if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'vcf'], sep: "\t")
    .map{ row-> tuple(row.name, file(row.vcf)) }
    .set { input_files }
}
else if (params.input_vcf) {
  input_vcf = file(params.input_vcf)
  Channel.fromList([tuple(input_vcf.name.take(input_vcf.name.lastIndexOf('.')), input_vcf)]).set { input_files }
}

workflow {
    if (params.filter) {
        FILTER_VCF(input_files)
        input_files = FILTER_VCF.out.filtered_vcfs
    }
    SUMMARY_VCF(input_files)
    NORMALIZE_VCF(input_files)
    SUMMARY_VCF_2(NORMALIZE_VCF.out.normalized_vcfs)
    NORMALIZE_VCF.out.normalized_vcfs
        .map {it.join("\t")}
        .collectFile(name: "${params.output}/normalized_vcfs.txt", newLine: true)
}

