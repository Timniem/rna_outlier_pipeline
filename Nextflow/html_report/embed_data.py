import json
import argparse
import pandas as pd
import plotly.express as px
import numpy as np
from jinja2 import Template
from gprofiler import GProfiler

#Plot functions
def scatter_plot_outrider(data, padj_treshold=0.05):
        """
        Scatter plot function for outrider (gene expression) data.
        --------------
        inppath_outpututs:
            data = pd.DataFrame

        outputs:
            fig_html = interactive plotly html code (HTML)    

        """
        data = data.replace([np.inf, -np.inf], 1)
        data = data.replace(np.nan, 1)
        data["minlogpVal"] = -np.log(data.pValue)
        data["significant"] = ['True' if padj < padj_treshold else "False" for padj in data["padjust"]]
        fig = px.scatter(data, x="zScore", y="minlogpVal", hover_data=["gene"],color_discrete_sequence=["rgba(255, 30, 30, 0.8)","rgba(60, 60, 60, 0.8)"], color='significant', labels={
                     "zScore": "zScore",
                     "minlogpVal": "-log pValue",
                 }, title="Expression volcano plot")
        fig_html = fig.to_html(full_html=False)
        return fig_html

def scatter_plot_fraser(data, padj_treshold=0.05):
        """
        Scatter plot function for outrider (gene expression) data.
        --------------
        inputs:
            data = pd.DataFrame

        outputs:
            fig_html = interactive plotly html code (HTML)    

        """
        data = data.replace([np.inf, -np.inf], 1)
        data = data.replace(np.nan, 1)
        data["minlogpVal"] = -np.log(data.pValue)
        data["significant"] = ['True' if padj < padj_treshold else "False" for padj in data["padjust"]]
        fig = px.scatter(data, x="deltaPsi", y="minlogpVal", hover_data=["gene"],color_discrete_sequence=["rgba(255, 30, 30, 0.8)","rgba(60, 60, 60, 0.8)"], color='significant', labels={
                     "deltaPsi": "deltaPsi",
                     "minlogpVal": "-log pValue",
                 }, title="Splicing volcano plot")
        fig_html = fig.to_html(full_html=False)
        return fig_html

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='RNA outlier HTML report module')
    parser.add_argument('-or','--outrider', help='Outrider results .tsv; example: path/to/outrider_results_patient._x.tsv', required=True)
    parser.add_argument('-fr','--fraser', help='Fraser results .tsv; example: path/to/fraser_results_patient_x.tsv', required=True)
    parser.add_argument('-t','--template', help='Html template; example: path/to/template.html', required=True)
    parser.add_argument('-s','--sampleid', help='Sample id used in reporting', required=True)
    parser.add_argument('-o','--output', help='Output path; example: path/to/patient_x_report.html', required=True)
    parser.add_argument('-g','--genepanel', help='Gene panel .txt file, one gene per line; example: path/to/genes.txt', required=False)
    parser.add_argument('-hp','--hpo', help='Human Phenotype Ontology .txt file, one HPO term per line; example: path/to/hpo.txt', required=False)
    args = vars(parser.parse_args())


    path_outr = args["outrider"] # example: path/to/outrider_results_patient._x.tsv
    path_frasr = args["fraser"] # example: path/to/fraser_results_patient_x.tsv
    html_template = args["template"] # example: path/to/template.html
    path_output = args["output"] # example: path/to/patient_x_report.html
    sample_id = args["sampleid"] # sampleId shown in report. 


    # load data in pandas
    outr_df = pd.read_csv(path_outr, sep="\t").rename(columns={"hgncSymbol":"gene"})[["gene","EnsemblID","pValue","padjust","zScore","l2fc", "rawcounts", "meanRawcounts", "normcounts", "meanCorrected"]]
    frasr_df = pd.read_csv(path_frasr, sep="\t").rename(columns={"hgncSymbol":"gene"})[["gene","chr", "start", "end", "width", "strand", "pValue","padjust","deltaPsi", "psiValue", "counts", "totalCounts"]]

    # create html plots
    outr_plot_html = scatter_plot_outrider(outr_df)
    frasr_plot_html = scatter_plot_fraser(frasr_df)

     # check for hpo and/or gene filters
    

    if args["genepanel"]:
           with open(args["genepanel"], 'r') as gene_file:
                genes = [line.strip().upper() for line in gene_file]
                outr_df = outr_df[outr_df["gene"].isin(genes)]
                frasr_df = frasr_df[frasr_df["gene"].isin(genes)]

    if args["hpo"]:
           with open(args["hpo"], 'r') as hpo_file:
                hpo_terms = [line.strip().upper() for line in hpo_file]
                hpo_df = pd.read_csv('resources/phenotype_to_genes.txt', sep='\t')
                hpo_df = hpo_df[["gene_symbol","hpo_id"]][hpo_df.hpo_id.isin(hpo_terms)]
                outr_df["hpo"] = [", ".join(hpo_df.hpo_id[hpo_df.gene_symbol == gene].unique().tolist()) for gene in outr_df["gene"]]
                frasr_df["hpo"] = [", ".join(hpo_df.hpo_id[hpo_df.gene_symbol == gene].unique().tolist()) for gene in frasr_df["gene"]]
    else:
        hpo_terms = False
        outr_df["hpo"] = ["" for gene in outr_df["gene"]]
        frasr_df["hpo"] = ["" for gene in frasr_df["gene"]]

    # reorder hpo to second position
    for df in [outr_df, frasr_df]:
        col = df.pop('hpo')
        df.insert(1, col.name, col)

    # hpo json object for active filters
    if hpo_terms:
          hpo_terms_json = json.dumps(hpo_terms)
    else:
          hpo_terms_json = json.dumps([])

    # create dataframe htmls
    outr_df_html = outr_df[outr_df["padjust"] < .99].to_html(index=False, classes='expression_table', border=0)
    frasr_df_html = frasr_df[frasr_df["padjust"] < .99].to_html(index=False, classes='splice_table', border=0)


    # gene enrichment analysis using gProfiler
    gp = GProfiler(return_dataframe=True)
    gene_query = [gene.split('.')[0] for gene in outr_df['EnsemblID'][outr_df["padjust"] < 0.05 ].tolist()]
    if gene_query:
        gene_enrichment_df = gp.profile(organism='hsapiens', query=gene_query)
        gene_enrichment_html = gene_enrichment_df.to_html(index=False, classes='gene_enrichment', border=0)
    else:
        gene_enrichment_html = " <h3>No result</h3>"

    # open html template
    with open(html_template, "r") as html_template_file:
           template = html_template_file.read()

    # Create a Jinja2 Template object
    jinja_template = Template(template)

    # Render the template with the table
    content = jinja_template.render({"expression_table": outr_df_html,
                                     "splice_table" : frasr_df_html,
                                     "plot_expression": outr_plot_html,
                                     "plot_splicing": frasr_plot_html,
                                     "gene_enrichment": gene_enrichment_html,
                                     "patient_id": sample_id,
                                     "hpo_terms": hpo_terms_json
                                     })

    # save the html report
    with open(path_output, 'w', encoding="utf-8") as output:
           output.write(content)
