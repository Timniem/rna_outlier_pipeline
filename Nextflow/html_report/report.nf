#!/usr/bin/env nextflow
/**
RNA outliers to interactive report workflow
**/

process ResultsToHtml {
    time '1h'
    memory '4 GB'
    cpus 1

    publishDir "$params.output/report/", mode: 'copy'

    input:
        tuple path(samplesheet), path(files)
    output:
        path "*.html"

    script:
        """
        ${CMD_REPORT} bash -c '
        which python3
        while IFS=\$"\t" read -r -a sampleid ; do
            if [ "\${sampleid}" != "sampleID" ]; then
                echo "Now processing \${sampleid}:"
                outrider_path="\${sampleid}_result_table_outrider.tsv"
                fraser_path="\${sampleid}_result_table_fraser.tsv"
                output_path="\${sampleid}_report.html"
                echo \$outrider_path
                echo \$fraser_path
                echo \$output_path
                # run html report script
                python3 ${params.report.embedScript} -or "\${outrider_path}" -fr "\${fraser_path}" -t "${params.report.htmlTemplate}" -s "\${sampleid}" -o "\${output_path}"
            fi
        done < $samplesheet
        '
        """
}