#!/bin/bash

output_folder=output/test0

nextflow main.nf --help

nextflow main.nf -profile test,conda --output $output_folder --skip_normalization

# missing SNpEff data folder
nextflow main.nf -profile test,conda --output $output_folder --snpeff_organism hg19
test ! -d $output_folder