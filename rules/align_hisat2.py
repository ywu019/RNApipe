import yaml 
import os
import pandas as pd
import numpy as np

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile, Loader=yaml.FullLoader)

if config["TRIM"] == "fastp":
    path = "/trim/fastp"
elif config["TRIM"] == "cutadapt":
    path = "/trim/cutadapt"

END=config["END"]
output_path = config["OUTPUTPATH"]
samples=pd.read_csv(config["SAMPLES"],dtype=str,sep='\t').set_index(["condition","replicate"],drop=False)
condition=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['condition']))
replicate=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['replicate']))
pairend=["R1","R2"]

rule all:
    input:
       sortbam=expand(output_path + "/align/hisat2/{condition}_{replicate}.sort.bam",condition=condition,replicate=replicate)

rule hisat2_index:
    input:
        genome=config["GENOME"],
        gtf=config["ANNOTATION"]
    output:
        dict_index=directory(os.path.dirname(config["GENOME"])+"/hisat2"),
        hisat2_splice=os.path.dirname(config["GENOME"])+"/hisat2/splicesites.txt"
    params:
        index=os.path.dirname(config["GENOME"])+"/hisat2/genome",
        thread=config["THREAD"]
    shell:
        """
        hisat2-build -p {params.thread} {input.genome} {params.index}
        hisat2_extract_splice_sites.py {input.gtf} > {output.hisat2_splice}
        """

if END=="pair":
    rule hisat2:
        input:
            TrimR1_fastq=lambda wildcards:expand(output_path+path+"/{condition}_{replicate}_R1_trimmed.fastq.gz",condition=wildcards.condition,replicate=wildcards.replicate),
            TrimR2_fastq=lambda wildcards:expand(output_path+path+"/{condition}_{replicate}_R2_trimmed.fastq.gz",condition=wildcards.condition,replicate=wildcards.replicate),
            #dict_index=directory(os.path.dirname(config["GENOME"])+"/hisat2"),
            hisat2_splice=os.path.dirname(config["GENOME"])+"/hisat2/splicesites.txt"
        output:
            sam=temp(output_path + "/align/hisat2/{condition}_{replicate}.sam"),
            bam=temp(output_path + "/align/hisat2/{condition}_{replicate}.bam")
        params:
            index=os.path.dirname(config["GENOME"])+"/hisat2/genome",
            summary=output_path + "/align/hisat2/{condition}_{replicate}.summary.txt",
            thread=config["THREAD"]
        shell:
            """
            hisat2 -p {params.thread} -x {params.index} -1 {input.TrimR1_fastq} -2 {input.TrimR2_fastq} -S {output.sam} --summary-file {params.summary}
            samtools view -@ {params.thread} -bS {output.sam} > {output.bam}
            """
else:
    rule hisat2:
        input:
            Trim_fastq=lambda wildcards:expand(output_path+path+"/{condition}_{replicate}_trimmed.fastq.gz",condition=wildcards.condition,replicate=wildcards.replicate),
            #dict_index=directory(os.path.dirname(config["GENOME"])+"/hisat2"),
            hisat2_splice=os.path.dirname(config["GENOME"])+"/hisat2/splicesites.txt"
        output:
            sam=temp(output_path + "/align/hisat2/{condition}_{replicate}.sam"),
            bam=temp(output_path + "/align/hisat2/{condition}_{replicate}.bam")
        params:
            index=os.path.dirname(config["GENOME"])+"/hisat2/genome",
            summary=output_path + "/align/hisat2/{condition}_{replicate}.summary.txt",
            thread=config["THREAD"]
        shell:
            """
            hisat2 -p {params.thread} -x {params.index} -U {input.Trim_fastq} -S {output.sam} --summary-file {params.summary}
            samtools view -@ {params.thread} -bS {output.sam} > {output.bam}
            """

rule sort:
    input:
        bam=output_path + "/align/hisat2/{condition}_{replicate}.bam"
    output:
        sortbam=output_path + "/align/hisat2/{condition}_{replicate}.sort.bam"
    params:
        thread=config["THREAD"]
    shell:
        """
        samtools sort -@ {params.thread} {input.bam} -o {output.sortbam}
        """
