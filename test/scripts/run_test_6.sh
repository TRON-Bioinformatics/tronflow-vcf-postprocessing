#!/bin/bash


source test/scripts/assert.sh
output_folder=test/output/test6


echo -e "tumor_normal\tprimary:"`pwd`"/test/data/TESTX_S1_L001.bam" > test/data/test_bams.txt
echo -e "tumor_normal\tnormal:"`pwd`"/test/data/TESTX_S1_L002.bam" >> test/data/test_bams.txt
echo -e "single_sample\ttumor:"`pwd`"/test/data/TESTX_S1_L001.bam" >> test/data/test_bams.txt
echo -e "single_sample\ttumor:"`pwd`"/test/data/TESTX_S1_L002.bam" >> test/data/test_bams.txt
echo -e "single_sample\tnormal:"`pwd`"/test/data/TESTX_S1_L001.bam" >> test/data/test_bams.txt
echo -e "single_sample\tnormal:"`pwd`"/test/data/TESTX_S1_L002.bam" >> test/data/test_bams.txt

# run
nextflow main.nf -profile test,conda,ci --output $output_folder --input_bams test/data/test_bams.txt

# check results
test -s $output_folder/single_sample/single_sample.normalized.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.normalized.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/single_sample/single_sample.vaf.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.vaf.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/single_sample/single_sample.filtered_multiallelics.vcf || { echo "Missing test 6 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.filtered_multiallelics.vcf || { echo "Missing test 6 output file!"; exit 1; }

assert_eq `grep -v '#' $output_folder/single_sample/single_sample.normalized.vcf | wc -l | cut -d' ' -f 1` 32 "Wrong number of variants"
assert_eq `grep -v '#' $output_folder/tumor_normal/tumor_normal.normalized.vcf | wc -l | cut -d' ' -f 1` 32 "Wrong number of variants"
assert_eq `grep -v '#' $output_folder/single_sample/single_sample.vaf.vcf | wc -l | cut -d' ' -f 1` 32 "Wrong number of variants"
assert_eq `grep tumor_af $output_folder/single_sample/single_sample.vaf.vcf | wc -l | cut -d' ' -f 1` 35 "Wrong number of variants"
assert_eq `grep normal_af $output_folder/single_sample/single_sample.vaf.vcf | wc -l | cut -d' ' -f 1` 35 "Wrong number of variants"
assert_eq `grep -v '#' $output_folder/tumor_normal/tumor_normal.vaf.vcf | wc -l | cut -d' ' -f 1` 32 "Wrong number of variants"
assert_eq `grep -v '#' $output_folder/tumor_normal/tumor_normal.filtered_multiallelics.vcf | wc -l | cut -d' ' -f 1` 23 "Wrong number of variants"
assert_eq `grep primary_af $output_folder/tumor_normal/tumor_normal.filtered_multiallelics.vcf | wc -l | cut -d' ' -f 1` 24 "Wrong number of variants"