/**
Nextflow FRASER workflow

**/
nextflow.enable.dsl=2


process FraserCount {
    time '8h'
    memory '16 GB'
    cpus 1

    input:
        path frasercountR
        path samplesheet
        path bamFiles
        path baiFiles
    output:
        path "fraser_output"

    script: 
        """
        mkdir "fraser_output"
        ${CMD_OUTRIDER_FRASER} Rscript "${frasercountR}" "${samplesheet}" "fraser_output"
        """
}


process MergeCounts {
    time '1h'
    memory '46 GB'
    cpus 1

    input:
        path ext_counts
        path mergescriptR
        path fraser_output
        val ext_amount_fraser
        

    output:

        path fraser_output

    script: 
        """
        ${CMD_OUTRIDER_FRASER} Rscript ${mergescriptR} "${fraser_output}" "${ext_counts}" "${ext_amount_fraser}"
        """

}

process Fraser {
    time '20h'
    memory '84 GB'
    cpus 6

    publishDir "$params.output/fraser", mode: 'copy'

    input:
        path samplesheet
        path output
        path fraserR

    output:

        path "*.tsv"

    script: 
        """
        ${CMD_OUTRIDER_FRASER} Rscript "${fraserR}" "${samplesheet}" "${output}"
        """
}