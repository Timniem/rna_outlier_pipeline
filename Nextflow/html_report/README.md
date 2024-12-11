<code>
  
  usage: embed_data.py [-h] -or OUTRIDER -fr FRASER -t TEMPLATE -s SAMPLEID -o OUTPUT [-g GENEPANEL] [-hp HPO]
  
  RNA outlier HTML report module
  
  options:
    -h, --help            show this help message and exit
    -or OUTRIDER, --outrider OUTRIDER
                          Outrider results .tsv; example: path/to/outrider_results_patient._x.tsv
    -fr FRASER, --fraser FRASER
                          Fraser results .tsv; example: path/to/fraser_results_patient_x.tsv
    -t TEMPLATE, --template TEMPLATE
                          Html template; example: path/to/template.html
    -s SAMPLEID, --sampleid SAMPLEID
                          Sample id used in reporting
    -o OUTPUT, --output OUTPUT
                          Output path; example: path/to/patient_x_report.html
    -g GENEPANEL, --genepanel GENEPANEL
                          Gene panel .txt file, one gene per line; example: path/to/genes.txt
    -hp HPO, --hpo HPO    Human Phenotype Ontology .txt file, one HPO term per line; example: path/to/hpo.txt

</code>
