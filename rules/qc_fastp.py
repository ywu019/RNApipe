import yaml 
import os
import pandas as pd
import numpy as np

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile, Loader=yaml.FullLoader)

adapter=config["ADAPTER"]
END=config["END"]
input_path = config["READSPATH"]
output_path = config["OUTPUTPATH"]
samples=pd.read_csv(config["SAMPLES"],dtype=str,sep='\t').set_index(["condition","replicate"],drop=False)
condition=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['condition']))
replicate=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['replicate']))
pairend=["R1","R2"]

rule all:
    input:
        expand(output_path + "/visualize/report/report_quality_control_after_trimming.html")

if END=="pair":
    rule fastp:
        input:
            R1_fastq = input_path + "/{condition}_{replicate}_R1.fastq.gz",
            R2_fastq = input_path + "/{condition}_{replicate}_R2.fastq.gz"
        output:
            TrimR1_fastq=output_path+"/trim/fastp/{condition}_{replicate}_R1_trimmed.fastq.gz",
            TrimR2_fastq=output_path+"/trim/fastp/{condition}_{replicate}_R2_trimmed.fastq.gz"
        params:
            adapter=adapter,
            thread=config["THREAD"]
        log:
            html=output_path+"/trim/fastp/{condition}_{replicate}.fastp.html",
            json=output_path+"/trim/fastp/{condition}_{replicate}.fastp.json"
        shell:
            """
            fastp -w {params.thread} -a {params.adapter} -i {input.R1_fastq} -I {input.R2_fastq} -o {output.TrimR1_fastq} -O {output.TrimR2_fastq} --length_required 20 -h {log.html} -j {log.json}
            """
    
    rule fastqc:
        input:
            Trim_fastq=output_path+"/trim/fastp/{condition}_{replicate}_{pairend}_trimmed.fastq.gz"
        output:
            html=output_path+"/trim/fastp/fastqc_after_trimming/{condition}_{replicate}_{pairend}_trimmed_fastqc.html"
        params:
            thread=config["THREAD"],
            outputpath=output_path+"/trim/fastp/fastqc_after_trimming"
        shell:
            "fastqc -o {params.outputpath} -t {params.thread}  {input.Trim_fastq}"
    
    rule summaryReport:
        input:
            html = expand(output_path+"/trim/fastp/fastqc_after_trimming/{condition}_{replicate}_{pairend}_trimmed_fastqc.html",condition=condition,replicate=replicate,pairend=pairend)
        output:
            report = output_path + "/visualize/report/report_quality_control_after_trimming.html"
        params:
            path = output_path + "/trim/fastp/fastqc_after_trimming"
        shell:
            "multiqc {params.path} --filename {output.report}"
else:
    rule fastp:
        input:
            fastq = input_path + "/{condition}_{replicate}.fastq.gz"
        output:
            Trim_fastq=output_path+"/trim/fastp/{condition}_{replicate}_trimmed.fastq.gz"
        params:
            adapter=adapter,
            thread=config["THREAD"]
        log:
            html=output_path+"/trim/fastp/{condition}_{replicate}.fastp.html",
            json=output_path+"/trim/fastp/{condition}_{replicate}.fastp.json"
        shell:
            """
            fastp -w {params.thread} -a {params.adapter} -i {input.fastq} -o {output.Trim_fastq} --length_required 20 -h {log.html} -j {log.json}
            """
    
    rule fastqc:
        input:
            Trim_fastq=output_path+"/trim/fastp/{condition}_{replicate}_trimmed.fastq.gz"
        output:
            html=output_path+"/trim/fastp/fastqc_after_trimming/{condition}_{replicate}_trimmed_fastqc.html"
        params:
            thread=config["THREAD"],
            outputpath=output_path+"/trim/fastp/fastqc_after_trimming"
        shell:
            "fastqc -o {params.outputpath} -t {params.thread}  {input.Trim_fastq}"
    
    rule summaryReport:
        input:
            html = expand(output_path+"/trim/fastp/fastqc_after_trimming/{condition}_{replicate}_trimmed_fastqc.html",condition=condition,replicate=replicate)
        output:
            report = output_path + "/visualize/report/report_quality_control_after_trimming.html"
        params:
            path = output_path + "/trim/fastp/fastqc_after_trimming"
        shell:
            "multiqc {params.path} --filename {output.report}"
    


