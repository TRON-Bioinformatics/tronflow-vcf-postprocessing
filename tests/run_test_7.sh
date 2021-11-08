#!/bin/bash


source bin/assert.sh
output_folder=output/test7
echo -e "tumor_normal\t"`pwd`"/test_data/TESTX_S1_L001.bam\t"`pwd`"/test_data/TESTX_S1_L002.bam" > test_data/test_bams.txt
echo -e "single_sample\t"`pwd`"/test_data/TESTX_S1_L001.bam,"`pwd`"/test_data/TESTX_S1_L002.bam\t"`pwd`"/test_data/TESTX_S1_L001.bam,"`pwd`"/test_data/TESTX_S1_L002.bam" >> test_data/test_bams.txt
nextflow main.nf -profile test,conda --output $output_folder --input_bams test_data/test_bams.txt --skip_multiallelic_filter
test -s $output_folder/single_sample/single_sample.normalized.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.normalized.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/single_sample/single_sample.normalized.vaf.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.normalized.vaf.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/single_sample/single_sample.normalized.vaf.filtered_multiallelics.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.normalized.vaf.filtered_multiallelics.vcf || { echo "Missing test 6 output file!"; exit 1; }
assert_eq `wc -l $output_folder/single_sample/single_sample.normalized.vcf | cut -d' ' -f 1` 53 "Wrong number of variants"
assert_eq `wc -l $output_folder/tumor_normal/tumor_normal.normalized.vcf | cut -d' ' -f 1` 53 "Wrong number of variants"
assert_eq `wc -l $output_folder/single_sample/single_sample.normalized.vaf.vcf | cut -d' ' -f 1` 72 "Wrong number of variants"
assert_eq `grep tumor_af $output_folder/single_sample/single_sample.normalized.vaf.vcf | wc -l | cut -d' ' -f 1` 35 "Wrong number of variants"
assert_eq `grep normal_af $output_folder/single_sample/single_sample.normalized.vaf.vcf | wc -l | cut -d' ' -f 1` 35 "Wrong number of variants"
assert_eq `wc -l $output_folder/tumor_normal/tumor_normal.normalized.vaf.vcf | cut -d' ' -f 1` 60 "Wrong number of variants"