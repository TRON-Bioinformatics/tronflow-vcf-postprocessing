#!/bin/bash

### this test cannot be automated as it relies on SnpEff references which need to be downloaded beforehand
source test/scripts/assert.sh
output_folder=test/output/test14
data_folder=`pwd`/test/data

nextflow main.nf -profile test,conda --output $output_folder --skip_normalization

test -s $output_folder/single_sample/single_sample.reheader.vcf || { echo "Missing test 14 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.reheader.vcf || { echo "Missing test 14 output file!"; exit 1; }


