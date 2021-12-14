import yaml 
import os
import pandas as pd
import numpy as np

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile)

if config["DEA"] == "DESeq2":
    path = "/DEA/DESeq2"
elif config["DEA"] == "edgeR":
    path = "/DEA/edgeR"
elif config["DEA"] == "limma":
    path = "/DEA/limma"

output_path = config["OUTPUTPATH"]
ncbi_gtf = config["ANNOTATION"]
control = config["CONTROL"][0]
treat = config["TREAT"][0]
speice_abbr = config["SPEICE_ABBREVIATION"]


rule all:
    input:
        kobas_out = output_path +"/annotation/go_kegg/kobas_deg_"+ control + "_" + treat + ".txt",
        All_EntrezID = output_path +"/annotation/go_kegg/GTF_EntrezID_geneid.txt",
        KEGG_plot=output_path +"/visualize/annotation/go_kegg/KEGG_"+ control + "_" + treat + ".pdf",
        GO_plot=output_path +"/visualize/annotation/go_kegg/GO_"+ control + "_" + treat + ".pdf",
        GO_barplot=output_path +"/visualize/annotation/go_kegg/GO_"+ control + "_" + treat + "_barplot.pdf",
        #GO_txt=output_path +"/annotation/go_kegg/GO_"+ control + "_" + treat + "_deg.txt",
        #KEGG_txt=output_path +"/annotation/go_kegg/KEGG_"+ control + "_" + treat + "_deg.txt",

rule deg_Symbol:
    input:
        deg = output_path + path + "/deg_" + control + "_" + treat + ".tsv"
    output:
        Gene_Symbol = temp(output_path +"/annotation/go_kegg/deg_geneid_"+ control + "_" + treat + ".txt")
    shell:
        "awk 'NR==1{{next}}{{print $1}}' {input.deg} > {output.Gene_Symbol}"

rule NCBIGTF_GetID:
    input:
        Gene_Symbol = output_path +"/annotation/go_kegg/deg_geneid_"+ control + "_" + treat + ".txt",
        ncbi_gtf = ncbi_gtf
    output:
        All_EntrezID = output_path +"/annotation/go_kegg/GTF_EntrezID_geneid.txt"
    shell:
        """
        sed '/^#/d' {input.ncbi_gtf} | awk -F '\t' '$3=="gene"{{print $NF}}' |  awk -F ';' '{{for (i=1;i<=NF;i++){{if ($i ~/gene_id|db_xref "GeneID/) {{print $i}}}}}}' | awk '{{print $2}}' | sed -n '{{N;s/\\n/\\t/p}}' | sed -e 's/\"//g' -e 's/;//g' -e 's/GeneID://g' > {output.All_EntrezID}
        """

rule deg_EntrezID:
    input:
        Gene_Symbol = output_path +"/annotation/go_kegg/deg_geneid_"+ control + "_" + treat + ".txt",
        All_EntrezID = output_path +"/annotation/go_kegg/GTF_EntrezID_geneid.txt"
    output:
        Gene_EntrezID = temp(output_path +"/annotation/go_kegg/deg_EntrezID_"+ control + "_" + treat + ".txt")
    shell:
        """
        awk 'NR==FNR{{a[$0]}}(NR>FNR && ($1 in a)){{print $2}}' {input.Gene_Symbol} {input.All_EntrezID} > {output.Gene_EntrezID}
        """

rule kobas:
    input:
        Gene_EntrezID = output_path +"/annotation/go_kegg/deg_EntrezID_"+ control + "_" + treat + ".txt"
    output:
        kobas_out = output_path +"/annotation/go_kegg/kobas_deg_"+ control + "_" + treat + ".txt" 
    params:
        speice_abbr = speice_abbr,
        kobas_path = os.getcwd()+"/software/kobas-3.0/"
    conda:
        "../envs/kobas_RNApipe.yaml"
    shell:
        """
        ./software/kobas-3.0/scripts/run_kobas.py  -i {input.Gene_EntrezID} -t id:ncbigene -s {params.speice_abbr} -k {params.kobas_path} -o {output.kobas_out}
        """

rule GOKEGG:
    input:
        kobas_out = output_path +"/annotation/go_kegg/kobas_deg_"+ control + "_" + treat + ".txt"
    output:
        #GO_txt=output_path +"/annotation/go_kegg/GO_"+ control + "_" + treat + "_deg.txt",
        #KEGG_txt=output_path +"/annotation/go_kegg/KEGG_"+ control + "_" + treat + "_deg.txt",
        KEGG_plot=output_path +"/visualize/annotation/go_kegg/KEGG_"+ control + "_" + treat + ".pdf",
        GO_plot=output_path +"/visualize/annotation/go_kegg/GO_"+ control + "_" + treat + ".pdf",
        GO_barplot=output_path +"/visualize/annotation/go_kegg/GO_"+ control + "_" + treat + "_barplot.pdf"
    shell:
        "Rscript scripts/GOKEGG.r {input.kobas_out} {output.KEGG_plot} {output.GO_plot} {output.GO_barplot}"
