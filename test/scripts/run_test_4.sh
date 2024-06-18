#!/bin/bash


source test/scripts/assert.sh
output_folder=test/output/test4
nextflow main.nf -profile test,conda,ci --input_vcfs test/data/test_input_no_ad.txt --output $output_folder
test -s $output_folder/sample_no_ad/sample_no_ad.normalized.vcf || { echo "Missing test 4 output file!"; exit 1; }
assert_eq `wc -l $output_folder/sample_no_ad/sample_no_ad.normalized.vcf | cut -d' ' -f 1` 52 "Wrong number of variants"