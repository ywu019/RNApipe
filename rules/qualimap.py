import yaml 
import os
import pandas as pd
import numpy as np

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile)

if config["ALIGN"] == "hisat2":
    path = "/align/hisat2"
elif config["ALIGN"] == "STAR":
    path = "/align/star"

if config["QUANTIFY"] == "featureCounts":
    count_path = "/quantify/featureCounts"
elif config["QUANTIFY"] == "htseq":
    count_path = "/quantify/htseq"

END=config["END"]
output_path = config["OUTPUTPATH"]
samples=pd.read_csv(config["SAMPLES"],dtype=str,sep='\t').set_index(["condition","replicate"],drop=False)
condition=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['condition']))
replicate=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['replicate']))

rule all:
    input:
        report =output_path + "/visualize/report/report_align_count.html"

rule alignmentQC:
    input:
        sortbam=output_path + path + "/{condition}_{replicate}.sort.bam"
    output:
        bamqc = temp(directory(output_path + path +"/{condition}_{replicate}_BAMqc"))
    shell:
        "qualimap bamqc -bam {input.sortbam} -nt 10 --java-mem-size=6G -outdir {output.bamqc}"


rule summaryReport:
    input:
        bamqc = expand(output_path + path +"/{condition}_{replicate}_BAMqc",condition=condition,replicate=replicate),
        count_summary = output_path + count_path + "/genome.tsv.summary"
    output:
        report =output_path + "/visualize/report/report_align_count.html"
    shell:
        "multiqc {input.bamqc} {input.count_summary} --filename {output.report}"
