#!/bin/bash


source test/scripts/assert.sh
output_folder=test/output/test12

# build input BAMs file
echo -e "tumor_normal\tprimary:"`pwd`"/test/data/TESTX_S1_L001.bam" > test/data/test_bams.txt
echo -e "tumor_normal\tnormal:"`pwd`"/test/data/TESTX_S1_L002.bam" >> test/data/test_bams.txt
echo -e "single_sample\ttumor:"`pwd`"/test/data/TESTX_S1_L001.bam" >> test/data/test_bams.txt
echo -e "single_sample\ttumor:"`pwd`"/test/data/TESTX_S1_L002.bam" >> test/data/test_bams.txt
echo -e "single_sample\tnormal:"`pwd`"/test/data/TESTX_S1_L001.bam" >> test/data/test_bams.txt
echo -e "single_sample\tnormal:"`pwd`"/test/data/TESTX_S1_L002.bam" >> test/data/test_bams.txt

# build input purities file
echo -e "tumor_normal\tprimary:0.5" > test/data/test_purities.txt
echo -e "tumor_normal\tnormal:0.6" >> test/data/test_purities.txt
echo -e "single_sample\ttumor:0.7" >> test/data/test_purities.txt
echo -e "single_sample\tnormal:0.8" >> test/data/test_purities.txt

# build input clonalities file
echo -e "tumor_normal\tprimary:`pwd`/test/data/best.seg.minimal.my_metastasis.bed" > test/data/test_clonalities.txt
echo -e "tumor_normal\tnormal:`pwd`/test/data/best.seg.minimal.my_tumor.bed" >> test/data/test_clonalities.txt
echo -e "single_sample\ttumor:3" >> test/data/test_clonalities.txt
echo -e "single_sample\tnormal:2" >> test/data/test_clonalities.txt

nextflow main.nf -profile test,conda --output $output_folder \
--input_bams test/data/test_bams.txt \
--input_purities test/data/test_purities.txt \
--input_clonalities test/data/test_clonalities.txt \
--skip_normalization

# test output files
test -s $output_folder/single_sample/single_sample.filtered_multiallelics.vcf || { echo "Missing test 10 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.filtered_multiallelics.vcf || { echo "Missing test 10 output file!"; exit 1; }
test -s $output_folder/single_sample/single_sample.vaf.vcf || { echo "Missing test 10 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.vaf.vcf || { echo "Missing test 10 output file!"; exit 1; }