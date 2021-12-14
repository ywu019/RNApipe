import yaml
import os
import pandas as pd
import numpy as np

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile)

if config["ALIGN"] == "hisat2":
    path = "/align/hisat2"
elif config["ALIGN"] == "STAR":
    path = "/align/STAR"

output_path = config["OUTPUTPATH"]
samples=pd.read_csv(config["SAMPLES"],dtype=str,sep='\t').set_index(["condition","replicate"],drop=False)
condition=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['condition']))
replicate=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['replicate']))

rule all:
    input:
        bigwig_file = expand(output_path + path + "/{condition}_{replicate}.bigwig",condition=condition,replicate=replicate),
        baifile = expand(output_path + path + "/{condition}_{replicate}.sort.bam.bai",condition=condition,replicate=replicate)

rule bam_index:
    input:
        sortbam=output_path + path + "/{condition}_{replicate}.sort.bam"
    output:
        baifile=output_path + path + "/{condition}_{replicate}.sort.bam.bai"
    shell:
        """
        samtools index {input.sortbam}
        """

rule bamtobigwig:
    input:
        sortbam=output_path + path + "/{condition}_{replicate}.sort.bam"
    output:
        bigwig_file=output_path + path + "/{condition}_{replicate}.bigwig"
    shell:
        """
        bamCoverage -p 10 -bs 1 --bam {input.sortbam} --outFileName {output.bigwig_file} --outFileFormat bigwig --normalizeUsing RPKM --exactScaling
        """


