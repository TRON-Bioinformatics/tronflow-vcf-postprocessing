vcf=$1


# keeps only the passed variants
bcftools view --apply-filters PASS -o `basename $vcf .vcf`.passed.vcf $vcf
