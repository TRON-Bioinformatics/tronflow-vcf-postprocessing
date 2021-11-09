#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_NORM; VT_DECOMPOSE_COMPLEX; REMOVE_DUPLICATES } from './modules/normalization'
include { FILTER_VCF } from './modules/filter'
include { SUMMARY_VCF; SUMMARY_VCF as SUMMARY_VCF_2 } from './modules/summary'
include { VAFATOR; MULTIALLELIC_FILTER } from './modules/vafator'

params.help= false
params.input_vcfs = false
params.input_bams = false
params.input_vcf = false
params.reference = false
params.output = "output"
params.skip_decompose_complex = false
params.filter = false
params.cpus = 1
params.memory = "4g"
params.vcf_without_ad = false
params.mapping_quality = false
params.base_call_quality = false
params.skip_multiallelic_filter = false
params.prefix = false


if (params.help) {
    log.info params.help_message
    exit 0
}

if (! params.input_vcfs && ! params.input_vcf) {
  exit 1, "Neither --input_vcfs or --input_vcf are provided!"
}
else if (params.input_vcfs && params.input_vcf) {
  exit 1, "Both --input_vcfs and --input_vcf are provided! Please, provide only one."
}
else if (params.input_vcfs) {
  Channel
    .fromPath(params.input_vcfs)
    .splitCsv(header: ['name', 'vcf'], sep: "\t")
    .map{ row-> tuple(row.name, file(row.vcf)) }
    .set{ input_vcfs }
}
else if (params.input_vcf) {
  input_vcf = file(params.input_vcf)
  Channel.fromList( [ tuple( input_vcf.name.take( input_vcf.name.lastIndexOf('.') ), input_vcf) ] ).set { input_vcfs }
}

if (params.input_bams) {
    Channel
    .fromPath(params.input_bams)
    .splitCsv(header: ['name', 'tumor_bams', 'normal_bams'], sep: "\t")
    .map{ row-> tuple(row.name, row.tumor_bams, row.normal_bams) }
    .set { input_bams }
}

workflow {

    if (params.filter) {
        FILTER_VCF(input_vcfs)
        input_vcfs = FILTER_VCF.out.filtered_vcfs
    }

    SUMMARY_VCF(input_vcfs)

    final_vcfs = BCFTOOLS_NORM(input_vcfs)
    if (! params.skip_decompose_complex) {
        VT_DECOMPOSE_COMPLEX(final_vcfs)
        final_vcfs = VT_DECOMPOSE_COMPLEX.out.decomposed_vcfs
    }
    REMOVE_DUPLICATES(final_vcfs)
    final_vcfs = REMOVE_DUPLICATES.out.deduplicated_vcfs

    SUMMARY_VCF_2(final_vcfs)

    if ( params.input_bams) {
        VAFATOR(final_vcfs.join(input_bams))
        final_vcfs = VAFATOR.out.annotated_vcf
        if ( ! params.skip_multiallelic_filter ) {
            final_vcfs = MULTIALLELIC_FILTER(final_vcfs)
            final_vcfs = MULTIALLELIC_FILTER.out.filtered_vcf
        }
    }

    final_vcfs.map {it.join("\t")}.collectFile(name: "${params.output}/normalized_vcfs.txt", newLine: true)
}

