import yaml 
import os
import pandas as pd
import numpy as np

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile, Loader=yaml.FullLoader)

if config["QUANTIFY"] == "featureCounts":
    path = "/quantify/featureCounts"
elif config["QUANTIFY"] == "htseq":
    path = "/quantify/htseq"

output_path = config["OUTPUTPATH"]
control = config["CONTROL"][0]
treat = config["TREAT"][0]
speice_abbr = config["SPEICE_ABBREVIATION"]
con = [treat, control]

rule all:
    input:
        #pre_txt = output_path + "/annotation/gsea/Gene_id_gsea.txt",
        #gct = output_path + "/annotation/gsea/gsea.txt",
        #cls = output_path + "/annotation/gsea/sample.cls",
        vis_pdf = expand(output_path + "/visualize/annotation/gsea/{con}.pdf", con=con)

rule get_pre_txt:
    input:
        tpm_tsv = output_path + path + "/genome_tpm.tsv"
    output:
        pre_txt = temp(output_path + "/annotation/gsea/Gene_id_gsea.txt")
    shell:
        """
        grep '0.0' -w -v {input.tpm_tsv} > {output.pre_txt}
        """

rule prepare_gsea_file:
    input:
        pre_txt = output_path + "/annotation/gsea/Gene_id_gsea.txt",
        EntrezID_list = output_path + "/annotation/go_kegg/GTF_EntrezID_geneid.txt",
        Samples = config["SAMPLES"]
    output:
        gct = temp(output_path + "/annotation/gsea/gsea.txt"),
        cls = temp(output_path + "/annotation/gsea/sample.cls")
    params:
        treat = treat,
        control = control
    shell:
        """
        python scripts/get_gsea_file.py -i {input.pre_txt} -s {input.Samples} --treat {params.treat} --control {params.control} -e {input.EntrezID_list} -g {output.gct} -c {output.cls}
        """

rule gsea:
    input:
        gct = output_path + "/annotation/gsea/gsea.txt",
        cls = output_path + "/annotation/gsea/sample.cls",
    output:
        pdf = expand(output_path + "/annotation/gsea/{con}.pdf", con=con)
    params:
        final_path = output_path + "/annotation/gsea/",
        species = config["SPEICE_ABBREVIATION"],
    conda:
        "../envs/kobas_RNApipe.yaml"
    shell:
        """
        ./software/kobas-3.0/scripts/run_gsea.py -e {input.gct} -c {input.cls} -i id:ncbigene -s {params.species} -k ./software/kobas-3.0/ -o {params.final_path}
        """

rule visualize:
    input:
        pdf = expand(output_path + "/annotation/gsea/{con}.pdf", con=con)
    output:
        vis_pdf = expand(output_path + "/visualize/annotation/gsea/{con}.pdf", con=con)
    params:
        final_path = output_path + "/visualize/annotation/gsea/"
    shell:
        """
        cp {input.pdf} {params.final_path}
        """
