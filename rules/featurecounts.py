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

END=config["END"]
output_path = config["OUTPUTPATH"]
samples=pd.read_csv(config["SAMPLES"],dtype=str,sep='\t').set_index(["condition","replicate"],drop=False)
condition=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['condition']))
replicate=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['replicate']))
pairend=["R1","R2"]

rule all:
    input:
        tpm_tsv=output_path + "/quantify/featureCounts/genome_tpm.tsv"

if END=='pair':
    rule featureCounts:
        input:
            sortbam=expand(output_path + path + "/{condition}_{replicate}.sort.bam",condition=condition,replicate=replicate),
            gtf=config["ANNOTATION"]
        output:
            count_tsv=output_path + "/quantify/featureCounts/genome.tsv"
        shell:
            """
            featureCounts -T 10 -p -t exon -g gene_id  -a {input.gtf} -o {output.count_tsv} {input.sortbam} &&\
            sed -i '1d' {output.count_tsv} 
            """
else:
    rule featureCounts:
        input:
            sortbam=expand(output_path + path + "/{condition}_{replicate}.sort.bam",condition=condition,replicate=replicate),
            gtf=config["ANNOTATION"]
        output:
            count_tsv=temp(output_path + "/quantify/featureCounts/genome.tsv")
        shell:
            """
            featureCounts -T 10 -t exon -g gene_id  -a {input.gtf} -o {output.count_tsv} {input.sortbam} &&\
            sed -i '1d' {output.count_tsv} 
            """

rule Counts2tpm:
    input:
        count_tsv=lambda wildcards:expand(output_path + "/quantify/featureCounts/genome.tsv")
    output:
        tpm_tsv=output_path + "/quantify/featureCounts/genome_tpm.tsv",
        count_tsv=output_path + "/quantify/featureCounts/genome_count.tsv"
    shell:
        "python scripts/featurecounts2tpm.py -i {input.count_tsv} -o1 {output.count_tsv} -o2 {output.tpm_tsv}"
