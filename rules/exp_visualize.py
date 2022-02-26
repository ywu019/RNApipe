import yaml 
import os

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile, Loader=yaml.FullLoader)

if config["QUANTIFY"] == "featureCounts":
    path = "/quantify/featureCounts"
elif config["QUANTIFY"] == "htseq":
    path = "/quantify/htseq"

control = config["CONTROL"][0]
treat = config["TREAT"][0]
output_path = config["OUTPUTPATH"]

rule all:
    input:
        plot_box = output_path + "/visualize/exp/violinplot_" + control + "_" + treat + ".pdf",
        plot_cumulative = output_path + "/visualize/exp/cumulative_" + control + "_" + treat + ".pdf",
        plot_pca = output_path + "/visualize/exp/PCA_" + control + "_" + treat + ".pdf",
        plot_correlation_heatmap = output_path + "/visualize/exp/correlation_" + control + "_" + treat + ".pdf",
        plot_fvid = output_path + "/visualize/exp/fvid_" + control + "_" + treat + ".pdf",
        plot_tpm_heatmap = output_path + "/visualize/exp/heatmap_" + control + "_" + treat + ".pdf"

rule visualize:
    input:
        tpm_tsv=output_path + path + "/genome_tpm.tsv"
    output:
        plot_box = output_path + "/visualize/exp/violinplot_" + control + "_" + treat + ".pdf",
        plot_cumulative = output_path + "/visualize/exp/cumulative_" + control + "_" + treat + ".pdf",
        plot_pca = output_path + "/visualize/exp/PCA_" + control + "_" + treat + ".pdf",
        plot_correlation_heatmap = output_path + "/visualize/exp/correlation_" + control + "_" + treat + ".pdf",
        plot_fvid = output_path + "/visualize/exp/fvid_" + control + "_" + treat + ".pdf",
        plot_tpm_heatmap = output_path + "/visualize/exp/heatmap_" + control + "_" + treat + ".pdf"
    shell:
        """
        Rscript scripts/exp_visualize.r {input.tpm_tsv} {output.plot_box} {output.plot_cumulative} {output.plot_pca} {output.plot_correlation_heatmap} {output.plot_fvid} {output.plot_tpm_heatmap}
        """
