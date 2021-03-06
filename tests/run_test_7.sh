#!/bin/bash


source tests/assert.sh
output_folder=output/test7
echo -e "tumor_normal\tprimary:"`pwd`"/test_data/TESTX_S1_L001.bam" > test_data/test_bams.txt
echo -e "tumor_normal\tnormal:"`pwd`"/test_data/TESTX_S1_L002.bam" >> test_data/test_bams.txt
echo -e "single_sample\ttumor:"`pwd`"/test_data/TESTX_S1_L001.bam" >> test_data/test_bams.txt
echo -e "single_sample\ttumor:"`pwd`"/test_data/TESTX_S1_L002.bam" >> test_data/test_bams.txt
echo -e "single_sample\tnormal:"`pwd`"/test_data/TESTX_S1_L001.bam" >> test_data/test_bams.txt
echo -e "single_sample\tnormal:"`pwd`"/test_data/TESTX_S1_L002.bam" >> test_data/test_bams.txt
nextflow main.nf -profile test,conda --output $output_folder --input_bams test_data/test_bams.txt --skip_multiallelic_filter
test -s $output_folder/single_sample/single_sample.normalized.vcf || { echo "Missing test 7 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.normalized.vcf || { echo "Missing test 7 output file!"; exit 1; }
test -s $output_folder/single_sample/single_sample.vaf.vcf || { echo "Missing test 7 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.vaf.vcf || { echo "Missing test 7 output file!"; exit 1; }
assert_eq `wc -l $output_folder/single_sample/single_sample.normalized.vcf | cut -d' ' -f 1` 53 "Wrong number of variants"
assert_eq `wc -l $output_folder/tumor_normal/tumor_normal.normalized.vcf | cut -d' ' -f 1` 53 "Wrong number of variants"
assert_eq `wc -l $output_folder/single_sample/single_sample.vaf.vcf | cut -d' ' -f 1` 72 "Wrong number of variants"
assert_eq `grep tumor_af $output_folder/single_sample/single_sample.vaf.vcf | wc -l | cut -d' ' -f 1` 35 "Wrong number of variants"
assert_eq `grep normal_af $output_folder/single_sample/single_sample.vaf.vcf | wc -l | cut -d' ' -f 1` 35 "Wrong number of variants"
assert_eq `wc -l $output_folder/tumor_normal/tumor_normal.vaf.vcf | cut -d' ' -f 1` 60 "Wrong number of variants"