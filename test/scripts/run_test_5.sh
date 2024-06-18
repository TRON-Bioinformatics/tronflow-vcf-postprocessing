#!/bin/bash


source test/scripts/assert.sh
output_folder=test/output/test5
nextflow main.nf -profile test,conda,ci --output $output_folder --input_vcfs false --input_vcf test/data/test_single_sample.vcf
test -s $output_folder/test_single_sample/test_single_sample.normalized.vcf || { echo "Missing test 4 output file!"; exit 1; }
assert_eq `wc -l $output_folder/test_single_sample/test_single_sample.normalized.vcf | cut -d' ' -f 1` 53 "Wrong number of variants"