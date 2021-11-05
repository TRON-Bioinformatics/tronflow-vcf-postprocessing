#!/bin/bash


source bin/assert.sh
output_folder=output/test2
nextflow main.nf -profile test,conda --output $output_folder --filter PASS,MNV
test -s $output_folder/single_sample/single_sample.normalized.vcf || { echo "Missing test 2 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.normalized.vcf || { echo "Missing test 2 output file!"; exit 1; }
assert `wc -l $output_folder/single_sample/single_sample.normalized.vcf` 35
assert `wc -l $output_folder/tumor_normal/tumor_normal.normalized.vcf` 35