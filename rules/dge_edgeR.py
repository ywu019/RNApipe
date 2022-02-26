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
conditions = np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['condition']))
control = config["CONTROL"][0]
treat = config["TREAT"][0]

rule all:
    input:
        volcano = output_path + "/visualize/DEA/edgeR/volcano_" + control + "_" + treat + ".pdf",
        heatmap = output_path + "/visualize/DEA/edgeR/heatmap_" + control + "_" + treat + ".pdf"

rule DEA:
    input:
        count_tsv=output_path + path + "/genome_count.tsv"
    output:
        dea = output_path + "/DEA/edgeR/dea_" + control + "_" + treat + ".tsv",
        deg = output_path + "/DEA/edgeR/deg_" + control + "_" + treat + ".tsv"
    shell:
        """
        Rscript scripts/dge_edgeR.r {input.count_tsv} {output.dea} {output.deg}
        sed 's/logFC/log2FoldChange/g' {output.dea} -i
        sed 's/logFC/log2FoldChange/g' {output.deg} -i
        sed 's/FDR/padj/g' {output.dea} -i
        sed 's/FDR/padj/g' {output.deg} -i
        """

rule visualize:
    input:
        dea = output_path + "/DEA/edgeR/dea_" + control + "_" + treat + ".tsv",
        tpm_tsv=output_path + path + "/genome_tpm.tsv"
    output:
        volcano = output_path + "/visualize/DEA/edgeR/volcano_" + control + "_" + treat + ".pdf",
        heatmap = output_path + "/visualize/DEA/edgeR/heatmap_" + control + "_" + treat + ".pdf"
    shell:
        "Rscript scripts/dge_visualize.r {input.dea} {input.tpm_tsv} {output.volcano} {output.heatmap}"
