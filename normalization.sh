#!/usr/bin/env bash

vcf=$1
reference=$2


module load anaconda/3/2019

# separate first by variant type
bcftools view --types snps -o `basename $vcf .vcf`.original.snps.vcf $vcf
bcftools view --types indels -o `basename $vcf .vcf`.original.indels.vcf $vcf
bcftools view --types mnps -o `basename $vcf .vcf`.original.mnps.vcf $vcf
bcftools view --types ref -o `basename $vcf .vcf`.original.ref.vcf $vcf
bcftools view --types bnd -o `basename $vcf .vcf`.original.bnd.vcf $vcf
bcftools view --types other -o `basename $vcf .vcf`.original.other.vcf $vcf

# decompose biallelic block substitutions (AC>TG to A>T and C>G)
# -a: best guess for non blocked substitutions
# -p: output phased genotypes and PS annotation
vt decompose_blocksub $vcf -a -p -o `basename $vcf .vcf`.atomic.vcf 2> `basename $vcf .vcf`.decompose_blocksub_stats.log

# decompose multiallelic variants into biallelic (C>T,G to C>T and C>G)
vt decompose `basename $vcf .vcf`.atomic.vcf -o `basename $vcf .vcf`.biallelic.vcf 2> `basename $vcf .vcf`.decompose_stats.log

# sort the input VCF
vt sort `basename $vcf .vcf`.biallelic.vcf -o `basename $vcf .vcf`.sorted.vcf

# normalize variants (trim and left alignment)
vt normalize `basename $vcf .vcf`.sorted.vcf -r $reference -o `basename $vcf .vcf`.normalized.vcf 2> `basename $vcf .vcf`.normalization_stats.log

# removes duplicated variants
vt uniq `basename $vcf .vcf`.normalized.vcf -o `basename $vcf .vcf`.uniq.vcf 2> `basename $vcf .vcf`.uniq_stats.log

# separate by variant type once normalized
bcftools view --types snps -o `basename $vcf .vcf`.normalized.snps.vcf `basename $vcf .vcf`.uniq.vcf
bcftools view --types indels -o `basename $vcf .vcf`.normalized.indels.vcf `basename $vcf .vcf`.uniq.vcf
bcftools view --types mnps -o `basename $vcf .vcf`.normalized.mnps.vcf `basename $vcf .vcf`.uniq.vcf
bcftools view --types ref -o `basename $vcf .vcf`.normalized.ref.vcf `basename $vcf .vcf`.uniq.vcf
bcftools view --types bnd -o `basename $vcf .vcf`.normalized.bnd.vcf `basename $vcf .vcf`.uniq.vcf
bcftools view --types other -o `basename $vcf .vcf`.normalized.other.vcf `basename $vcf .vcf`.uniq.vcf

# delete intermediate files
rm -f `basename $vcf .vcf`.atomic.vcf
rm -f `basename $vcf .vcf`.biallelic.vcf
rm -f `basename $vcf .vcf`.sorted.vcf
rm -f `basename $vcf .vcf`.normalized.vcf
rm -f `basename $vcf .vcf`.uniq.vcf
