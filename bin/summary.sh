#!/usr/bin/env bash
vcf=$1

module load anaconda/3/2019

bcftools summary $vcf > `basename $vcf .vcf`.stats
plot-vcfstats -p `basename $vcf .vcf`.stats_plots `basename $vcf .vcf`.stats
