#!/bin/bash


source test/scripts/assert.sh
output_folder=test/output/test10

echo -e "tumor_normal\tprimary:"`pwd`"/test/data/TESTX_S1_L001.bam" > test/data/test_bams.txt
echo -e "tumor_normal\tnormal:"`pwd`"/test/data/TESTX_S1_L002.bam" >> test/data/test_bams.txt
echo -e "single_sample\ttumor:"`pwd`"/test/data/TESTX_S1_L001.bam" >> test/data/test_bams.txt
echo -e "single_sample\ttumor:"`pwd`"/test/data/TESTX_S1_L002.bam" >> test/data/test_bams.txt
echo -e "single_sample\tnormal:"`pwd`"/test/data/TESTX_S1_L001.bam" >> test/data/test_bams.txt
echo -e "single_sample\tnormal:"`pwd`"/test/data/TESTX_S1_L002.bam" >> test/data/test_bams.txt

nextflow main.nf -profile test,conda --output $output_folder --input_bams test/data/test_bams.txt --skip_normalization

test -s $output_folder/single_sample/single_sample.filtered_multiallelics.vcf || { echo "Missing test 10 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.filtered_multiallelics.vcf || { echo "Missing test 10 output file!"; exit 1; }
test -s $output_folder/single_sample/single_sample.vaf.vcf || { echo "Missing test 10 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.vaf.vcf || { echo "Missing test 10 output file!"; exit 1; }