#!/bin/bash

### this test cannot be automated as it relies on SnpEff references which need to be downloaded beforehand
source test/scripts/assert.sh
output_folder=test/output/test13
data_folder=`pwd`/test/data

nextflow main.nf -profile test,conda --output $output_folder \
  --gff ${data_folder}/Homo_sapiens.GRCh37.87.minimal.gff3.gz \
  --reference ${data_folder}/ucsc.hg19.minimal.fasta

test -s $output_folder/single_sample/single_sample.annotated.vcf || { echo "Missing test 1 output file!"; exit 1; }
test -s $output_folder/tumor_normal/tumor_normal.annotated.vcf || { echo "Missing test 1 output file!"; exit 1; }

assert_eq `grep -v '#' $output_folder/single_sample/single_sample.annotated.vcf| grep BCSQ | wc -l | cut -d' ' -f 1` 32 "Wrong number of variants"
assert_eq `grep -v '#' $output_folder/tumor_normal/tumor_normal.annotated.vcf | grep BCSQ | wc -l | cut -d' ' -f 1` 32 "Wrong number of variants"