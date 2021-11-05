#!/bin/bash


source bin/assert.sh
output_folder=output/test3
nextflow main.nf -profile test,conda --output $output_folder --skip_decompose_complex
test -s $output_folder/single_sample/single_sample.normalized.vcf || { echo "Missing test 2 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.normalized.vcf || { echo "Missing test 2 output file!"; exit 1; }
assert `wc -l $output_folder/single_sample/single_sample.normalized.vcf` 48
assert `wc -l $output_folder/tumor_normal/tumor_normal.normalized.vcf` 48