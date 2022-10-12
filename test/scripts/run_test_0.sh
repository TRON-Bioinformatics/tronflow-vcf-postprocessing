#!/bin/bash

output_folder=test/output/test0

nextflow main.nf --help

# missing SNpEff data folder
nextflow main.nf -profile test,conda --output $output_folder --snpeff_organism hg19
test ! -d $output_folder