import yaml 
import os
import pandas as pd
import numpy as np
#os.chdir('/home/zpzhou/wuyue/snakemake/RNApipe')
with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile)

END=config["END"]
input_path = config["READSPATH"]
intermediate_path = config["OUTPUTPATH"]
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] 

samples=pd.read_csv(config["SAMPLES"],dtype=str,sep='\t').set_index(["condition","replicate"],drop=False)
condition=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['condition']))
replicate=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['replicate']))
pairend=["R1","R2"]


rule all:
    input:
        expand(intermediate_path + "/fastqc/report_quality_control.html")

if END=="pair":
    rule fastqc:
        input:
            fastq=expand(input_path + "/{condition}_{replicate}_{pairend}.fastq.gz",condition=condition,replicate=replicate,pairend=pairend)
        output:
            html=temp(intermediate_path+"/fastqc/{condition}_{replicate}_{pairend}_fastqc.html")
        params:
            threads=10,
            outputpath=intermediate_path+"/fastqc"
        shell:
            "fastqc -o {params.outputpath} -t {params.threads}  {input.fastq}"

    rule summaryReport:
        input:
            fastq = expand(intermediate_path+"/fastqc/{condition}_{replicate}_{pairend}_fastqc.html",condition=condition,replicate=replicate,pairend=pairend)
        output:
            report = intermediate_path + "/fastqc/report_quality_control.html"
        params:
            path = intermediate_path + "/fastqc"
        shell:
            "multiqc {params.path} --filename {output.report}"
else:
    rule fastqc:
        input:
            fastq=expand(input_path + "/{condition}_{replicate}.fastq.gz",condition=condition,replicate=replicate)
        output:
            html=temp(intermediate_path+"/fastqc/{condition}_{replicate}_fastqc.html")
        params:
            threads=10,
            outputpath=intermediate_path+"/fastqc"
        shell:
            "fastqc -o {params.outputpath} -t {params.threads}  {input.fastq}"

    rule summaryReport:
        input:
            fastq = expand(intermediate_path+"/fastqc/{condition}_{replicate}_fastqc.html",condition=condition,replicate=replicate)
        output:
            report = intermediate_path + "/fastqc/report_quality_control.html"
        params:
            path = intermediate_path + "/fastqc"
        shell:
            "multiqc {params.path} --filename {output.report}"

