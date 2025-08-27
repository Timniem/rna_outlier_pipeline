import sys
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
    fig = px.scatter(data, x="zScore", y="minlogpVal", hover_data=["gene"], color_discrete_map={"True":"rgba(255, 30, 30, 0.8)","False":"rgba(60, 60, 60, 0.8)"}, color='significant', labels={
                    "zScore": "zScore",
                    "minlogpVal": "-log pValue",
                }, title="Expression volcano plot")
    fig.update_layout(
                        plot_bgcolor='rgba(0, 0, 0, 0)',
                        paper_bgcolor='rgba(0, 0, 0, 0)',
                        width=500,
                        title_x=0.5
                    ),
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
    fig = px.scatter(data, x="deltaPsi", y="minlogpVal", hover_data=["gene"], color_discrete_map={"True":"rgba(255, 30, 30, 0.8)","False":"rgba(60, 60, 60, 0.8)"}, color='significant', labels={
                    "deltaPsi": "deltaPsi",
                    "minlogpVal": "-log pValue",
                }, title="Splicing volcano plot")
    fig.update_layout(
                        plot_bgcolor='rgba(0, 0, 0, 0)',
                        paper_bgcolor='rgba(0, 0, 0, 0)',
                        width=500,
                        title_x=0.5
                    ),
    fig_html = fig.to_html(full_html=False)
    return fig_html

def scatter_plot_mae(data, padj_treshold=0.05):
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
    fig = px.scatter(data, x="log2FC", y="minlogpVal", hover_data=["gene"], color_discrete_map={"True":"rgba(255, 30, 30, 0.8)","False":"rgba(60, 60, 60, 0.8)"}, color='significant', labels={
                    "log2FC": "log2FC",
                    "minlogpVal": "-log pValue",
                }, title="MAE volcano plot")
    fig.update_layout(
                        plot_bgcolor='rgba(0, 0, 0, 0)',
                        paper_bgcolor='rgba(0, 0, 0, 0)',
                        width=500,
                        title_x=0.5
                    ),
    fig_html = fig.to_html(full_html=False)
    return fig_html

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='RNA outlier HTML report module')
    parser.add_argument('-or','--outrider', help='Outrider results .tsv; example: path/to/outrider_results_patient._x.tsv', required=True)
    parser.add_argument('-fr','--fraser', help='Fraser results .tsv; example: path/to/fraser_results_patient_x.tsv', required=True)
    parser.add_argument('-ma','--mae', help='MAE results .tsv; example: path/to/mae_results_patient_x.tsv', required=False)
    parser.add_argument('-t','--template', help='Html template; example: path/to/template.html', required=True)
    parser.add_argument('-s','--sampleid', help='Sample id used in reporting', required=True)
    parser.add_argument('-o','--output', help='Output path; example: path/to/patient_x_report.html', required=True)
    parser.add_argument('-g','--genepanel', help='Gene panel .txt file, one gene per line; example: path/to/genes.txt', required=False)
    parser.add_argument('-hp','--hpo', help='Human Phenotype Ontology .txt file, one HPO term per line; example: path/to/hpo.txt', required=False)
    parser.add_argument('-hpf','--hpoFile', help='Human Phenotype Ontology phenotype_to_genes file example: path/to/phenotype_to_genes.txt', required=False)
    args = vars(parser.parse_args())


    path_outr = args["outrider"] # example: path/to/outrider_results_patient._x.tsv
    path_frasr = args["fraser"] # example: path/to/fraser_results_patient_x.tsv
    path_mae = args["mae"] # example: path/to/fraser_results_patient_x.tsv
    html_template = args["template"] # example: path/to/template.html
    path_output = args["output"] # example: path/to/patient_x_report.html
    sample_id = args["sampleid"] # sampleId shown in report.

    if not ".tsv" in args["mae"]: #MAE can be empty
        args["mae"] = None

    # load data in pandas
    outr_df = pd.read_csv(path_outr, sep="\t").rename(columns={"hgncSymbol":"gene"})[["gene","EnsemblID","pValue","padjust","zScore","l2fc", "rawcounts", "meanRawcounts", "normcounts", "meanCorrected"]]
    frasr_df = pd.read_csv(path_frasr, sep="\t").rename(columns={"hgncSymbol":"gene"})[["gene","chr", "start", "end", "width", "strand", "pValue","padjust","deltaPsi", "psiValue", "counts", "totalCounts"]]
    if args["mae"]:
        mae_df = pd.read_csv(path_mae, sep="\t").rename(columns={"hgncSymbol":"gene", "pvalue":"pValue", "end":"pos"})[["gene", "chr", "pos","refAllele" ,"altAllele" ,"pValue", "padjust", "log2FC","refCount","altCount","totalCount"]]
        mae_plot_html = scatter_plot_mae(mae_df)
    else:
        mae_plot_html = "plot not available"
    # create html plots
    outr_plot_html = scatter_plot_outrider(outr_df)
    frasr_plot_html = scatter_plot_fraser(frasr_df)
    

     # check for hpo and/or gene filters
    
    if args["genepanel"]:
        gene_panel_name = args["genepanel"]
        with open(args["genepanel"], 'r') as gene_file:
            genepanel_list = [line.strip().upper() for line in gene_file]
    else:
        gene_panel_name = 'none'
        genepanel_list = []

    if args["hpo"]:
        if args["hpoFile"]:
            path_hpo_file = args["hpoFile"] # example: path/to/phenotype_to_genes.txt
            with open(args["hpo"], 'r') as hpo_file:
                hpo_terms = [line.strip().upper() for line in hpo_file]
                hpo_df = pd.read_csv(path_hpo_file, sep='\t')
                hpo_df = hpo_df[["gene_symbol","hpo_id"]][hpo_df.hpo_id.isin(hpo_terms)]
                outr_df["hpo"] = [", ".join(hpo_df.hpo_id[hpo_df.gene_symbol == gene].unique().tolist()) for gene in outr_df["gene"]]
                frasr_df["hpo"] = [", ".join(hpo_df.hpo_id[hpo_df.gene_symbol == gene].unique().tolist()) for gene in frasr_df["gene"]]
                if args["mae"]:
                    mae_df["hpo"] = [", ".join(hpo_df.hpo_id[hpo_df.gene_symbol == gene].unique().tolist()) for gene in mae_df["gene"]]
        else:
            print("Please provide Human Phenotype Ontology phenotype_to_genes.txt file arg= -hpf --hpoFile")
            raise FileNotFoundError
    else:
        hpo_terms = False
        
        outr_df["hpo"] = ["" for gene in outr_df["gene"]]
        frasr_df["hpo"] = ["" for gene in frasr_df["gene"]]
        if args["mae"]:
            mae_df["hpo"] = ["" for gene in mae_df["gene"]]

    # reorder hpo to second position
    for df in [outr_df, frasr_df]:
        col = df.pop('hpo')
        df.insert(1, col.name, col)
    
    if args["mae"]:
        col = mae_df.pop('hpo')
        mae_df.insert(1, col.name, col)

    # hpo json object for active filters
    if hpo_terms:
        hpo_terms_json = json.dumps(hpo_terms)
    else:
        hpo_terms_json = json.dumps([])
    if genepanel_list:
        gene_list_json = json.dumps(genepanel_list)
    else:
        gene_list_json = json.dumps([])

    # gene enrichment analysis using gProfiler
    gp = GProfiler(return_dataframe=True)
    gene_query = [gene.split('.')[0] for gene in outr_df['EnsemblID'][outr_df["padjust"] < 0.05 ].tolist()]
    if gene_query:
        gene_enrichment_df = gp.profile(organism='hsapiens', query=gene_query)
        gene_enrichment_html = gene_enrichment_df.to_html(index=False, classes='gene_enrichment', border=0)
    else:
        gene_enrichment_html = "No results"

    # open html template
    with open(html_template, "r") as html_template_file:
           template = html_template_file.read()
    
    # create dataframe htmls
    if genepanel_list:
        outr_df = outr_df[outr_df["gene"].isin(genepanel_list)]
        frasr_d = frasr_df[frasr_df["gene"].isin(genepanel_list)]
        if args["mae"]:
            mae_df = mae_df[mae_df["gene"].isin(genepanel_list)]

    # create dataframe htmls
    outr_df_html = outr_df[outr_df["padjust"] < .99].to_html(index=False, classes='expression_table', border=0)
    frasr_df_html = frasr_df[frasr_df["padjust"] < .99].to_html(index=False, classes='splice_table', border=0)
    if args["mae"]:
        mae_df_html = mae_df[mae_df["padjust"] < .99].to_html(index=False, classes='mae_table', border=0)
    else:
        mae_df_html = "Table not available"

    # Create a Jinja2 Template object
    jinja_template = Template(template)

    # Render the template with the table
    content = jinja_template.render({"expression_table": outr_df_html,
                                     "splice_table" : frasr_df_html,
                                     "mae_table": mae_df_html,
                                     "expression_plot": outr_plot_html,
                                     "splice_plot": frasr_plot_html,
                                     "mae_plot": mae_plot_html,
                                     "gene_enrichment": gene_enrichment_html,
                                     "patient_id": sample_id,
                                     "hpo_terms": hpo_terms_json,
                                     "gene_panel_file": gene_panel_name,
                                     "R_version" : "4.3.1",
                                     "outrider_version": "1.20.1",
                                     "fraser_version": "1.99.4",
                                     "GATK_version": "4.2.3.0",
                                     "bcftools_version": "1.19",
                                     "tMAE_version": "mumichae/tMAE@1.0.0",
                                     "python_version": sys.version

                                     })

    # save the html report
    with open(path_output, 'w', encoding="utf-8") as output:
           output.write(content)
