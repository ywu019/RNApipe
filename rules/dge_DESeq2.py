import yaml 
import os
import pandas as pd
import numpy as np

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile)

if config["QUANTIFY"] == "featureCounts":
    path = "/quantify/featureCounts"
elif config["QUANTIFY"] == "htseq":
    path = "/quantify/htseq"

output_path = config["OUTPUTPATH"]
conditions = np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['condition']))
control = config["CONTROL"][0]
treat = config["TREAT"][0]

rule all:
    input:
        volcano = output_path + "/visualize/DEA/DESeq2/volcano_" + control + "_" + treat + ".pdf",
        heatmap = output_path + "/visualize/DEA/DESeq2/heatmap_" + control + "_" + treat + ".pdf"

rule DEA:
    input:
        count_tsv=output_path + path + "/genome_count.tsv"
    output:
        dea = output_path + "/DEA/DESeq2/dea_" + control + "_" + treat + ".tsv",
        deg = output_path + "/DEA/DESeq2/deg_" + control + "_" + treat + ".tsv"
    shell:
        "Rscript scripts/dge_DESeq2.r {input.count_tsv} {output.dea} {output.deg}"

rule visualize:
    input:
        dea = output_path + "/DEA/DESeq2/dea_" + control + "_" + treat + ".tsv",
        tpm_tsv=output_path + path + "/genome_tpm.tsv"
    output:
        volcano = output_path + "/visualize/DEA/DESeq2/volcano_" + control + "_" + treat + ".pdf",
        heatmap = output_path + "/visualize/DEA/DESeq2/heatmap_" + control + "_" + treat + ".pdf"
    shell:
        "Rscript scripts/dge_visualize.r {input.dea} {input.tpm_tsv} {output.volcano} {output.heatmap}"

