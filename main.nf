#!/usr/bin/env nextflow

params.help= false
params.input_files = false
params.reference = false
params.output = false
params.skip_split_mnps = false
params.filter = false
params.decompose_non_blocked_substitutions = false
params.skip_duplication_removal = false
params.skip_split_vcf_by_type = false
params.cpus = 1
params.memory = "4g"


if (params.help) {
    log.info params.help
    exit 0
}

publish_dir = "output"
if (params.output) {
  publish_dir = params.output
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'vcf'], sep: "\t")
    .map{ row-> tuple(row.name, file(row.vcf)) }
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}

if (params.filter) {
  process filterVcf {
    cpus params.cpus
    memory params.memory
    tag "${name}"


    input:
    	set name, file(vcf) from input_files

    output:
      set name, file("${vcf.baseName}.filtered.vcf") into filtered_vcf

    """
    # filter variants
    bcftools view --apply-filter ${params.filter} -o ${vcf.baseName}.filtered.vcf ${vcf}
    """
  }
}
else {
  filtered_vcf = input_files
}

/*
This step sets MAPQ to 0 for all unmapped reads + avoids soft clipping beyond the end of the reference genome
This step reorders chromosomes in the BAM file according to the provided reference (this step is required for GATK)
Adds the required read groups fields to the BAM file. The provided type is added to the BAM sample name.
*/
process normalizeVcf {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${publish_dir}/${name}", mode: "copy"

    input:
    	set name, file(vcf) from filtered_vcf

    output:
      set name, file("${name}.normalized.vcf") into normalized_vcf2
      set name, val("${publish_dir}/${name}/${name}.normalized.vcf") into normalized_vcf
      file("${name}.normalized.snps.vcf") optional true into normalized_snps_vcf_file
      file("${name}.normalized.indels.vcf") optional true into normalized_indels_vcf_file
      file("${name}.normalized.mnps.vcf") optional true into normalized_mnps_vcf_file
      file("${name}.normalized.bnd.vcf") optional true into normalized_bnd_vcf_file
      file("${name}.normalized.other.vcf") optional true into normalized_other_vcf_file

    script:
    decomposeNonBlockedSubstitutionsOption = params.decompose_non_blocked_substitutions ? " -a " : ""
    sortedVcf00 =  name + ".00.sorted.vcf"
    atomicVcf01 =  name + ".01.atomic.vcf"
    biallelicVcf02 =  name + ".02.biallelic.vcf"
    sortedVcf03 =  name + ".03.sorted.vcf"
    leftAlignedVcf04 =  name + ".04.left_aligned.vcf"
    normalizedVcf =  name + ".normalized.vcf"
    normalizedSnpsVcf =  name + ".normalized.snps.vcf"
    normalizedIndelsVcf =  name + ".normalized.indels.vcf"
    normalizedMnvsVcf =  name + ".normalized.mnps.vcf"
    normalizedOtherVcf =  name + ".normalized.other.vcf"
    normalizedBndVcf =  name + ".normalized.bnd.vcf"
    """
    # initial sort of the VCF
    vt sort ${vcf} -o ${sortedVcf00}

    # decompose biallelic block substitutions (AC>TG to A>T and C>G)
    # -a: best guess for non blocked substitutions
    # -p: output phased genotypes and PS annotation
    if (${params.skip_split_mnps}) ; then
      cp $sortedVcf00 ${atomicVcf01}
    else
      vt decompose_blocksub ${sortedVcf00} ${decomposeNonBlockedSubstitutionsOption} -o ${atomicVcf01}
    fi

    # decompose multiallelic variants into biallelic (C>T,G to C>T and C>G)
    vt decompose ${atomicVcf01} -o ${biallelicVcf02}

    # sort the input VCF
    vt sort ${biallelicVcf02} -o ${sortedVcf03}

    # normalize variants (trim and left alignment)
    vt normalize ${sortedVcf03} -r ${params.reference} -o ${leftAlignedVcf04}

    # removes duplicated variants
    if (${params.skip_duplication_removal}) ; then
        cp ${leftAlignedVcf04} ${normalizedVcf}
    else
        vt uniq ${leftAlignedVcf04} -o ${normalizedVcf}
    fi
    
    # separate by variant type once normalized
    if ( ! ${params.skip_split_vcf_by_type}) ; then
        # excludes other than SNP to avoid removing somatic reference variants
        bcftools view --exclude-types indels,mnps,bnd,other -o ${normalizedSnpsVcf} ${normalizedVcf}
        bcftools view --types indels -o ${normalizedIndelsVcf} ${normalizedVcf}
        bcftools view --types mnps -o ${normalizedMnvsVcf} ${normalizedVcf}
        bcftools view --types bnd -o ${normalizedBndVcf} ${normalizedVcf}
        bcftools view --types other -o ${normalizedOtherVcf} ${normalizedVcf}  
    fi

    # delete intermediate files
    rm -f ${atomicVcf01}
    rm -f ${biallelicVcf02}
    rm -f ${sortedVcf03}
    rm -f ${leftAlignedVcf04}
    """
}

process summaryVcf {
  cpus params.cpus
  memory params.memory
  tag "${name}"
  publishDir "${publish_dir}/${name}", mode: "copy"

  input:
    set name, file(vcf) from normalized_vcf2

  output:
    file("${name}_stats/*") into vcf_stats_plots

  """
  mkdir -p ${name}_stats
  bcftools stats $vcf > ${name}_stats/${vcf.baseName}.stats
  #plot-vcfstats -p ${name}_stats --no-PDF --title ${name} ${name}_stats/${vcf.baseName}.stats
  """
}

normalized_vcf
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/normalized_vcfs.txt", newLine: true)
