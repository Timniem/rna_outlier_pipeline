#!/usr/bin/env nextflow
/**
Nextflow MAE workflow
Uses Gagneur-lab MAE derived implementation from the DROP pipeline.. 
Results will be stored in a resultsfile (*.tsv)
**/

nextflow.enable.dsl=2


process MAEreadCounting {
    time '12h'
    memory '8 GB'
    cpus 1

    publishDir "$params.output/MAE/ASEreadcounts", mode: 'copy'

    input:
        tuple val(sampleID), path(vcf), path(bamFile)
        tuple val(fasta), path(fastafolder)

    output:
        val "${sampleID}"
        path "${sampleID}_maecounts.tsv"

    script:
    """
    ${CMD_MAE} bash -c '
    set -e
    set -u

    # From DROP:
    if samtools view -H "${bamFile}" | grep -q "@RG"; then
        printf "BAM contains read groups RG, continuing with ASEReadCounter...\n"
    else
        printf "%s\n" "ERROR: BAM file doesnt contain Read Group Tag RG"
        printf "RG doesnt exist, it can be added using:\n"
        printf "   gatk AddOrReplaceGroups -R /path/to/reference -I /your/input.bam -O /your/output.bam --QUIET true\n"
        printf " https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-\n"
        printf " Try rerunning this module using the BAM with RG tags\n"
        exit 1
    fi

    # Remove INFO field, it sometimes contains whitespace or invalid line lengths
    bcftools annotate -x "INFO" "$vcf" -o ${sampleID}_noinfo.vcf

    # Create index
    gatk IndexFeatureFile -I "${sampleID}_noinfo.vcf"

    # Only select heterozygous variants and only SNPs
    gatk SelectVariants -V ${sampleID}_noinfo.vcf -O "${sampleID}_temp.vcf" -R "${fasta}" --restrict-alleles-to BIALLELIC -select "vc.getHetCount()==1" --select-type-to-include SNP

    # Remove duplicates
    bcftools norm --rm-dup all "${sampleID}_temp.vcf" > "${sampleID}_temp.vcf.tmp" && mv "${sampleID}_temp.vcf.tmp" "${sampleID}_temp.vcf"

    # Create new index
    gatk IndexFeatureFile -I "${sampleID}_temp.vcf"

    # Run ASEReadCounter
    gatk ASEReadCounter -I "${bamFile}" -V "${sampleID}_temp.vcf" -O "${sampleID}_temp_counts.tsv" -R "${fasta}"

    # Sorting the count file
    (head -n 1 "${sampleID}_temp_counts.tsv" && tail -n +2 "${sampleID}_temp_counts.tsv" | sort -k1,1 -V -s) > "${sampleID}_maecounts.tsv"
    '
    """
}

process GetMAEresults {

    time '1h'
    memory '4 GB'
    cpus 1

    publishDir "$params.output/MAE/DeSEQresults", mode: 'copy'

    input:
        val sampleid
        path asecounts
        path resultsR

    output:
        
        path "${sampleid}_result_mae.tsv"

    script:
        """
        ${CMD_MAE} Rscript "${resultsR}" "${asecounts}" "${sampleid}_result_mae.tsv" "${sampleid}"
        """

}


