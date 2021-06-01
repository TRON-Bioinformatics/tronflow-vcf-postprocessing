# TronFlow variant normalization pipeline

[![DOI](https://zenodo.org/badge/372133189.svg)](https://zenodo.org/badge/latestdoi/372133189)

This pipeline aims at normalizing variants represented in a VCF into the convened normal form as described in Tan 2015. 
The variant normalization is based on the implementation in vt (Tan 2015) and bcftools (Danecek 2021). 
The pipeline is implemented on the Nextflow (Di Tommaso 2017) framework.
 
The pipeline consists of the following steps:
 * Variant filtering (optional)
 * Decomposition of MNPs into atomic variants (ie: AC > TG is decomposed into two variants A>T and C>G) (optional).
 * Decomposition of multiallelic variants into biallelic variants (ie: A > C,G is decomposed into two variants A > C and A > G)
 * Trim redundant sequence and left align indels, indels in repetitive sequences can have multiple representations
 * Remove duplicated variants
 
The output consists of:
 * The normalized VCF
 * Summary statistics before and after normalization


![Pipeline](images/variant_normalization_pipeline.png)

## Examples

### Variants are trimmed removing redundant bases

Before:
```
chr1	13083	.	AGA	AGC	.	UNTRIMMED		GT:AD	0:276,0	0/1:260,18
```
After:
```
chr1	13085	.	A	C	.	UNTRIMMED	OLD_CLUMPED=chr1|13083|AGA|AGC	GT:AD	0:276,276	0/1:260,260
```

**NOTE**: when a variant is changed during normalization the original variant is kept in the `INFO` field `OLD_CLUMPED`.

### Indels are left aligned

If an indel lies within a repetitive sequence it can be reported on different positions, the convention is to report 
the left most indel. This is what is called left alignment.

Before:
```
chr1	13141	.	C	AAAAAC	.	UNALIGNED		GT:AD	0:80,1	0/1:76,14
chr1	13141	.	CTGAGG	G	.	UNALIGNED		GT:AD	0:80,1	0/1:76,14
```

After:
```
chr1	13140	.	C	CAAAAA	.	UNALIGNED	OLD_CLUMPED=chr1|13141|C|AAAAAC	GT:AD	0:80,1	0/1:76,14
chr1	13140	.	CCTGAG	C	.	UNALIGNED	OLD_CLUMPED=chr1|13141|CTGAGG|G	GT:AD	0:80,1	0/1:76,14
```

**NOTE**: the rule of thumb to confirm if any given indel is left aligned, assuming that they are trimmed, is that the last base in both reference and 
alternate must differ.

### Multi-allelic variants are split.

Before:
```
chr1	13204	.	C	G,T	.	MULTIALLELIC		GT:AD	1:229,1,1	0/1:196,24,1 
```
After:
```
chr1	13204	.	C	G	.	MULTIALLELIC	OLD_VARIANT=chr1|13204|C|G,|1	GT:AD	1:229,229	0/1:196,196
chr1	13204	.	C	T	.	MULTIALLELIC	OLD_VARIANT=chr1|13204|C|G,|2	GT:AD	0:229,229	0/0:196,196
```

Note that the AD values are incorrectly set after the split, this seems to be a regression issue in bcftools v1.12 
and has been reported here https://github.com/samtools/bcftools/issues/1499.

### Multi Nucleotide Variants (MNVs) can be decomposed

Before:
```
chr1	13261	.	GCTCCT	CCCCCC	.	MNV		GT:AD	0:41,0	0/1:35,3
```

After:
```
chr1	13261	.	G	C	.	MNV	OLD_CLUMPED=chr1:13261:GCTCCT/CCCCCC	GT:AD:PS	0:41,0:13261	0|1:35,3:13261
chr1	13263	.	T	C	.	MNV	OLD_CLUMPED=chr1:13261:GCTCCT/CCCCCC	GT:AD:PS	0:41,0:13261	0|1:35,3:13261
chr1	13266	.	T	C	.	MNV	OLD_CLUMPED=chr1:13261:GCTCCT/CCCCCC	GT:AD:PS	0:41,0:13261	0|1:35,3:13261
```

This behaviour is optional and can be disabled with `--skip_decompose_complex`. Beware, that the phase is maintained in
the fields `FORMAT/GT` and `FORMAT/PS` as described in the VCF specification section 1.4.2.

### Complex variants combining SNV and indels can be decomposed

Before:
```
chr1	13321	.	AGCCCT	CGCC	.	MNV-INDEL		GT:AD	0:229,1	0/1:196,24
```

After:
```
chr1	13321	.	A	C	.	MNV-INDEL	OLD_CLUMPED=chr1:13321:AGCCCT/CGCC	GT:AD:PS	0:229,1:13321	0|1:196,24:13321
chr1	13322	.	GC	G	.	MNV-INDEL	OLD_CLUMPED=chr1:13321:AGCCCT/CGCC	GT:AD:PS	0:229,1:13321	0|1:196,24:13321
chr1	13325	.	CT	C	.	MNV-INDEL	OLD_CLUMPED=chr1:13321:AGCCCT/CGCC	GT:AD:PS	0:229,1:13321	0|1:196,24:13321
```

Same as MNVs this behaviour can de disabled with `--skip_decompose_complex`.



## How to run it

 ```
 $ nextflow run tron-bioinformatics/tronflow-variant-normalization --help
 
 TronFlow VCF normalization v${VERSION}

Usage:
    nextflow run main.nf --input_files input_files --reference reference.fasta


Input:
    * --input_vcf: the path to a single VCF to normalize (not compatible with --input_files)
    * --input_files: the path to a tab-separated values file containing in each row the sample name  and path to the VCF file (not compatible with --input_vcf)
    The input file does not have header!
    Example input file:
    sample1	/path/to/your/file.vcf
    sample2	/path/to/your/file2.vcf
    * --reference: path to the FASTA genome reference (indexes expected *.fai, *.dict)
    * --vcf-without-ad: indicate when the VCFs to normalize do not have the FORMAT/AD annotation

Optional input:
    * --output: the folder where to publish output
    * --skip_decompose_complex: flag indicating not to split complex variants (ie: MNVs and combinations of SNVs and indels)
    * --filter: specify the filter to apply if any (e.g.: PASS), only variants with this value will be kept

Output:
    * Normalized VCF file
    * Tab-separated values file with the absolute paths to the preprocessed BAMs, preprocessed_bams.txt
    * Summary statistics before and after normalization
 ```
 

## References

* Adrian Tan, Gonçalo R. Abecasis and Hyun Min Kang. Unified Representation of Genetic Variants. Bioinformatics (2015) 31(13): 2202-2204](http://bioinformatics.oxfordjournals.org/content/31/13/2202) and uses bcftools [Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics (Oxford, England), 27(21), 2987–2993. 10.1093/bioinformatics/btr509
* Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.
* Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. 10.1038/nbt.3820
