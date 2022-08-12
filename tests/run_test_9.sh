#!/bin/bash

### this test cannot be automated as it relies on SnpEff references which need to be downloaded beforehand
source tests/assert.sh
output_folder=output/test8
snpeff_datadir=/home/you/snpeff

nextflow main.nf -profile test,conda --output $output_folder --snpeff_organism hg19 --snpeff_datadir $snpeff_datadir

test -s $output_folder/single_sample/single_sample.normalized.vcf || { echo "Missing test 1 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.normalized.vcf || { echo "Missing test 1 output file!"; exit 1; }

assert_eq `grep -v '#' $output_folder/single_sample/single_sample.normalized.vcf | wc -l | cut -d' ' -f 1` 32 "Wrong number of variants"
assert_eq `grep -v '#' $output_folder/tumor_normal/tumor_normal.normalized.vcf | wc -l | cut -d' ' -f 1` 32 "Wrong number of variants"