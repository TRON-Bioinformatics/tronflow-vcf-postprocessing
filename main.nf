#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FILTER_VCF } from './modules/01_filter'
include { BCFTOOLS_NORM; VT_DECOMPOSE_COMPLEX; REMOVE_DUPLICATES } from './modules/02_normalization'
include { SUMMARY_VCF; SUMMARY_VCF as SUMMARY_VCF_2 } from './modules/03_summary'
include { VAFATOR; MULTIALLELIC_FILTER } from './modules/04_vafator'
include { VARIANT_ANNOTATION } from './modules/05_variant_annotation'

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
params.snpeff_organism = false
params.snpeff_args = ""
params.snpeff_datadir = false


if (params.help) {
    log.info params.help_message
    exit 0
}

if ( params.snpeff_organism && ! params.snpeff_datadir) {
  exit 1, "To run snpEff, please, provide your snpEff data folder with --snpeff_datadir"
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
    .splitCsv(header: ['name', 'sample_name', 'bam'], sep: "\t")
    .map{ row-> tuple(row.name, row.sample_name, row.bam) }
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

    if ( params.input_bams ) {
        VAFATOR(final_vcfs.join(input_bams.groupTuple()))
        final_vcfs = VAFATOR.out.annotated_vcf
        if ( ! params.skip_multiallelic_filter ) {
            final_vcfs = MULTIALLELIC_FILTER(final_vcfs)
            final_vcfs = MULTIALLELIC_FILTER.out.filtered_vcf
        }
    }

    if (params.snpeff_organism) {
        VARIANT_ANNOTATION(final_vcfs)
        final_vcfs = VARIANT_ANNOTATION.out.annotated_vcf
    }
}

